function SD = waveSourcePoints(cube, gridMask, time_ms, varargin)
% waveSourcePoints   Source/evaluation-point wave detection on an LFP cube.
%
% PURPOSE
% -------
% Branch-3 of the traveling-wave pipeline, ported from the T32
% wave_detection_sourcepoints.m and generalized over the trigger. Localizes wave
% ORIGINS on the grid frame-by-frame and scores them with a principled statistic,
% rather than a single global coherence. Per epoch, within a detection window:
%   1. bandpass 5-40 Hz -> generalized_phase (gp_lp=0, matches Muller demo)
%   2. phase-gradient vector field (dx,dy)
%   3. evaluation points  = moments the spatial-mean phase crosses `evaluation_angle`
%   4. source point       = argmax divergence of (dx,dy) at each ep
%   5. rho                = phase_correlation_distance(phase, source)
%   6. detected           = rho > rhoThr
%
% Disc-array masking (our 64-electrode octagon in an 8x10 frame): holes are filled
% per frame for gradient continuity, but the fabricated corners are kept OUT of
% detection -- NaN'd in the phase cube for the eval-point spatial mean (find_eval_
% points nansum keeps the angle correct) and for rho (phase_correlation_distance
% strips NaNs), and the source argmax is restricted to real electrodes.
%
% INPUTS
% ------
% cube      nEpoch x nRow x nCol x T (NaN at non-electrode cells).
% gridMask  nRow x nCol logical (true = real electrode).
% time_ms   1 x T epoch time axis (ms; 0 = trigger), from accumulateSaccadeCube.
% Name/value:
%   'Fs' 1000, 'filt_band' [5 40], 'filt_order' 4, 'gp_lp' 0,
%   'detect_win_ms' [-100 200], 'evaluation_angle' pi, 'ep_tol' 0.2,
%   'rhoThr' 0.5, 'pixelSpacing_mm' 0.7, 'Verbose' true
%   'NullDraws' 0     collect this many null rho samples by shuffling phase values
%                     across real electrodes at the evaluation points (Davis/Muller's
%                     null). 0 = off. Returned in SD.nullRho so a threshold can be
%                     CALIBRATED rather than assumed -- the default rhoThr 0.5 was
%                     never justified against a null on this array, and the paper's
%                     0.3 was calibrated on Utah-array noise, not ECoG.
%   'NullSeed'  5     rng seed for reproducible nulls
%   'SourceMask' []   nRow x nCol logical restricting where a source may be
%                     localized (default = gridMask; pass an eroded/interior mask
%                     to exclude rim electrodes as a boundary-artifact control)
%   'ReturnPhase' false  also return SD.phase{s} = nRow x nCol x nEP, the corner-
%                     masked generalized phase at each evaluation point. This is the
%                     EXACT frame rho was computed on, handed back rather than
%                     recomputed. The motif analysis clusters these frames, and a
%                     second bandpass+GP pass outside this function is a second bug
%                     surface -- waveAnalyticFrames documents four renderers that
%                     drifted from the detector that way. Off by default: it costs
%                     ~nEP*nRow*nCol doubles per epoch (a few MB per session).
%
% OUTPUT (struct SD)
% ------------------
%   .nDetections    nEpoch x 1
%   .maxRho         nEpoch x 1
%   .firstDetect_ms nEpoch x 1
%   .source_mean_rc nEpoch x 2  (rho-weighted mean source, row col)
%   .phase          nEpoch x 1 cell of nRow x nCol x nEP phase maps, ONLY when
%                   'ReturnPhase' is true (empty cell array otherwise)
%   .ep_idx .source_rc .rho .detected .dirDeg   nEpoch x 1 cell (per ep)
%                   dirDeg = propagation direction, top-CW (0=up, CW+), as waveGP.
%                   Send through waveDirToScreen before plotting.
%   .nullRho        pooled null rho samples when NullDraws > 0
%   .Params
%
% Requires wave-matlab + generalized-phase (+ smoothn in third_party) + Circular
% Statistics Toolbox on path.
%
% See also: waveGP, accumulateSaccadeCube

p = struct('Fs',1000, 'filt_band',[5 40], 'filt_order',4, 'gp_lp',0, ...
           'detect_win_ms',[-100 200], 'evaluation_angle',pi, 'ep_tol',0.2, ...
           'rhoThr',0.5, 'pixelSpacing_mm',0.7, 'SourceMask',[], ...
           'NullDraws',0, 'NullSeed',5, 'ReturnPhase',false, 'Verbose',true);
for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end
for fn = {'find_evaluation_points','find_source_points','smoothn','phase_correlation_distance'}
    if isempty(which(fn{1})), error('waveSourcePoints:dep','Missing on path: %s', fn{1}); end
end

checkPitchDefault('waveSourcePoints', varargin, gridMask);
[nE, nRow, nCol, T] = size(cube);
gmRC = logical(gridMask);
srcMask = gmRC; if ~isempty(p.SourceMask), srcMask = logical(p.SourceMask); end  % source-candidate cells
[Xg, Yg] = meshgrid(1:nCol, 1:nRow);
detMask = time_ms >= p.detect_win_ms(1) & time_ms <= p.detect_win_ms(2);
detIdx  = find(detMask); detIdxOffset = detIdx(1) - 1;

nDetections    = zeros(nE,1);
maxRho         = nan(nE,1);
firstDetect_ms = nan(nE,1);
source_mean_rc = nan(nE,2);
[ep_idx,source_rc,rho_c,detected_c,dir_c] = deal(cell(nE,1));
phase_c = cell(nE,1);
nullRho = [];
if p.NullDraws > 0, rng(p.NullSeed); end

for s = 1:nE
    raw = squeeze(cube(s,:,:,:));
    for k = 1:T                                    % fill holes (continuity only)
        fr = raw(:,:,k); if any(isnan(fr(:))), raw(:,:,k) = fillDiscFrame(fr); end
    end
    xf  = bandpass_filter(double(raw), p.filt_band(1), p.filt_band(2), p.filt_order, p.Fs);
    xgp = generalized_phase(xf, p.Fs, p.gp_lp);
    [~,~,dx,dy] = phase_gradient_complex_multiplication(xgp, p.pixelSpacing_mm);

    % masked phase cube: NaN the fabricated corners (per frame)
    xgpM = xgp;
    for k = 1:T, fr = xgpM(:,:,k); fr(~gmRC) = NaN; xgpM(:,:,k) = fr; end

    % evaluation points within the detection window
    ep_local = find_evaluation_points(xgpM(:,:,detMask), p.evaluation_angle, p.ep_tol);
    if isempty(ep_local)
        ep_idx{s}=[]; source_rc{s}=zeros(0,2); rho_c{s}=[]; detected_c{s}=false(0,1); dir_c{s}=[];
        phase_c{s}=zeros(nRow,nCol,0); continue;
    end
    ep_abs = ep_local(:) + detIdxOffset;

    % source points (masked divergence argmax -> real electrodes only)
    src = masked_source_points(ep_abs, Xg, Yg, dx, dy, srcMask); % 2 x nEP (x=col, y=row)

    % ---- per-event PROPAGATION DIRECTION ----
    % Uses propagationField, the single source of truth for the sign-corrected
    % propagation vectors (waveGP uses it too), so event directions cannot drift
    % from the rest of the suite. Amplitude-weighted spatial mean at the evaluation
    % point, then the same math->topCW conversion waveGP applies; downstream display
    % must still go through waveDirToScreen.
    Uprop = propagationField(xgp, gmRC, p.Fs);
    Amp   = abs(xgp);

    % rho per ep (NaN corners auto-stripped inside phase_correlation_distance)
    nEP = numel(ep_abs); rho = nan(1,nEP); epDirDeg = nan(1,nEP);
    for jj = 1:nEP
        if any(isnan(src(:,jj))), continue; end
        rho(jj) = phase_correlation_distance(angle(xgpM(:,:,ep_abs(jj))), src(:,jj), p.pixelSpacing_mm);

        u = Uprop(:,:,ep_abs(jj)); a = Amp(:,:,ep_abs(jj));
        u(~gmRC) = NaN; a(~gmRC) = NaN;
        g = sum(u(:).*a(:),'omitnan') / (sum(a(:),'omitnan') + eps);
        epDirDeg(jj) = wrap180deg(90 - rad2deg(angle(g)));      % top-CW, as waveGP
    end

    % Hand back the SAME corner-masked phase frames rho was just computed on. Indexing
    % xgpM at ep_abs reproduces angle(xgpM(:,:,ep_abs(jj))) exactly, element for element,
    % so a motif clustered on SD.phase is clustered on the detector's own view.
    if p.ReturnPhase, phase_c{s} = angle(xgpM(:,:,ep_abs)); end

    % ---- NULL for rho: shuffle phase values across real electrodes ----
    % Davis/Muller build the rho null by "randomly shuffling phase values on the
    % electrode array". Done HERE, inside the same loop, so the null passes through
    % the identical source/rho code path as the observed value -- a null computed
    % by a separate reimplementation is not a null, it is a second bug surface.
    % The source point is held fixed and only the phase field is permuted, which is
    % exactly the "is phase organised about THIS point beyond chance" question.
    if p.NullDraws > 0 && nEP > 0
        nPerEp = max(1, ceil(p.NullDraws / max(nE,1) / max(nEP,1)));
        for jj = 1:nEP
            if any(isnan(src(:,jj))), continue; end
            fr = angle(xgpM(:,:,ep_abs(jj)));
            ok = find(gmRC & isfinite(fr));
            if numel(ok) < 8, continue; end
            for q = 1:nPerEp
                frS = fr;
                frS(ok) = fr(ok(randperm(numel(ok))));
                nullRho(end+1) = phase_correlation_distance(frS, src(:,jj), p.pixelSpacing_mm); %#ok<AGROW>
            end
        end
    end

    detected = rho > p.rhoThr;
    nDetections(s) = sum(detected);
    if any(detected)
        rd = rho(detected); maxRho(s) = max(rd);
        firstDetect_ms(s) = time_ms(ep_abs(find(detected,1,'first')));
        w = rd(:)./sum(rd); srcD = src(:,detected);        % (col, row)
        source_mean_rc(s,:) = [srcD(2,:)*w, srcD(1,:)*w];   % (row, col)
    end
    ep_idx{s} = ep_abs; source_rc{s} = [src(2,:); src(1,:)].'; rho_c{s} = rho(:); detected_c{s} = detected(:); dir_c{s} = epDirDeg(:);
    if p.Verbose && mod(s,50)==0, fprintf('  waveSourcePoints %d/%d\n', s, nE); end
end

SD = struct('nDetections',nDetections, 'maxRho',maxRho, 'firstDetect_ms',firstDetect_ms, ...
            'source_mean_rc',source_mean_rc, 'ep_idx',{ep_idx}, 'source_rc',{source_rc}, ...
            'rho',{rho_c}, 'detected',{detected_c}, 'dirDeg',{dir_c}, 'phase',{phase_c}, ...
            'nullRho',nullRho(:), 'Params',p);

if p.Verbose && ~isempty(nullRho)
    fprintf(['  waveSourcePoints null: %d draws | rho null p50 %.3f p95 %.3f p99 %.3f ' ...
             '(current rhoThr %.2f)\n'], numel(nullRho), prctile(nullRho,50), ...
             prctile(nullRho,95), prctile(nullRho,99), p.rhoThr);
end
if p.Verbose
    fprintf('waveSourcePoints: %d epochs, fracDetected=%.2f, median(maxRho|det)=%.2f\n', ...
        nE, mean(nDetections>0), median(maxRho(nDetections>0),'omitnan'));
end
end

% ======================================================================
function src = masked_source_points(ep_abs, Xg, Yg, dx, dy, gmRC)
% find_source_points with the argmax restricted to real electrodes.
ws = warning('off','all'); cleaner = onCleanup(@() warning(ws)); %#ok<NASGU> % quiet smoothn/divergence
d = nan(size(Xg,1), size(Xg,2), numel(ep_abs));
for ii = 1:numel(ep_abs)
    d(:,:,ii) = divergence(Xg, Yg, dx(:,:,ep_abs(ii)), dy(:,:,ep_abs(ii)));
end
d = smoothn(d, 'robust');
src = nan(2, numel(ep_abs));
for ii = 1:numel(ep_abs)
    di = d(:,:,ii); di(~gmRC) = -Inf;              % no source on a fabricated cell
    [yy,xx] = find(di == max(di(:)));
    if numel(yy) == 1, src(:,ii) = [xx; yy]; end   % (x=col, y=row)
end
end

function a = wrap180deg(a)
a = mod(a + 180, 360) - 180;
end
