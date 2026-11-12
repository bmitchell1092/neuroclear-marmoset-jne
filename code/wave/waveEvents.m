function [E, S, P] = waveEvents(cube, gridMask, time_ms, varargin)
% waveEvents   Statistically-thresholded discrete traveling-wave events.
%
% PURPOSE
% -------
% Turns an LFP cube into a table of WAVE EVENTS: discrete moments at which the
% phase field is organised about a source point beyond what shuffled phase
% produces. Implements Davis, Muller et al. (Nature 2020) via waveSourcePoints,
% with the detection threshold CALIBRATED on this array's own null rather than
% inherited.
%
% WHY NOT R_dir (waveGP)
% ----------------------
% R_dir is an amplitude-weighted spatial coherence and rises with SNR. Measured on
% this array (2026-08-17), its post-saccade rise survives a spatial shuffle almost
% intact -- only ~4% excess -- so on its own it is weak evidence for propagation.
% The circular-linear correlation of phase with distance from a source point (rho)
% behaves completely differently: observed p50 0.33-0.44 against a null p50 of
% 0.15, with observed p95 exceeding the null p99 by 0.27-0.39. Use rho for
% detection; R_dir remains useful as a descriptive time course, not as evidence.
%
% THRESHOLD
% ---------
% 'RhoThr' 'auto' calibrates from the null: threshold = the 'NullPct' percentile
% (default 99) of shuffled-phase rho. On two light-rest sessions this landed at
% 0.368 and 0.369 -- reassuringly stable, and near the 0.3 that Davis et al.
% calibrated on Utah-array noise despite our ECoG having ~4x their inter-channel
% correlation (0.86-0.95 vs 0.2-0.4). Pass a number to override.
%
% INPUTS
% ------
% cube      nEpoch x nRow x nCol x T (NaN at non-electrode cells)
% gridMask  nRow x nCol logical
% time_ms   1 x T (ms, 0 = trigger)
% Name/value:
%   'RhoThr'    'auto'   'auto' | numeric
%   'NullPct'   99       percentile of the null used when RhoThr='auto'
%   'NullDraws' 4000     null samples to collect for the calibration
%   'InteriorOnly' true  restrict sources to interior electrodes (rim electrodes
%                        produce boundary artifacts in the divergence argmax)
%   'Verbose'   true
%   ... any waveSourcePoints option (detect_win_ms, filt_band, ...) is forwarded.
%       detect_win_ms defaults HERE to the full span of time_ms. waveSourcePoints'
%       own default is [-100 200], which silently discards the rest of the epoch --
%       see the note at the forwarding site.
%
% OUTPUTS
% -------
% E  event table, one row per detection: epoch, time_ms, rho, srcRow, srcCol, dirDeg
%    dirDeg is the propagation direction in the pipeline's native top-CW frame
%    (0 = up, CW+). Pass it through waveDirToScreen before any polar plot.
% S  summary struct: .rhoThr .nullPct .nullRho .nEpochs .nEval .nEvents
%    .detectRate (fraction of evaluation points that are events)
%    .eventsPerSec (per epoch-second inside the detection window)
%    .epochsWithEvent (fraction of epochs carrying >=1 event)
% P  nRow x nCol x height(E) phase maps (rad, NaN off gridMask), ONE PER EVENT ROW and
%    in the same order -- P(:,:,k) is the spatial phase pattern of event E(k,:). Empty
%    unless this third output is requested; asking for it switches 'ReturnPhase' on in
%    waveSourcePoints, so the map is the detector's own frame rather than a recompute.
%    This is the input the motif analysis clusters: an EVENT and a motif PATTERN are the
%    same object seen twice, which is what lets a motif time course be laid directly over
%    the event-rate time course without a method mismatch.
%
% See also: waveSourcePoints, waveGP, accumulateSaccadeCube, waveMotifCluster

p = struct('RhoThr','auto', 'NullPct',99, 'NullDraws',4000, ...
           'InteriorOnly',true, 'Verbose',true);
known = fieldnames(p); fwd = {};
for i = 1:2:numel(varargin)
    if any(strcmpi(varargin{i}, known)), p.(varargin{i}) = varargin{i+1};
    else, fwd(end+1:end+2) = varargin(i:i+1); end
end

% DEFAULT THE DETECTION WINDOW TO THE WHOLE EPOCH.
% waveSourcePoints defaults detect_win_ms to [-100 200] -- a stimulus-locked window
% baked into a general engine. Events simply DO NOT EXIST outside it, silently, so a
% caller that passes a 700 ms cube and never thinks about detect_win_ms gets events
% over 300 ms of it and rate statistics computed against a mostly-empty baseline.
% That is exactly what happened to run_motionWaveEvents (2026-09-02): its pre-window
% [-250 -50] sat almost entirely outside the default range, so the "pre" rate was
% divided by a window most of which could not contain a detection, and the driver
% reported a spurious 3.4x evoked rate ratio. Detecting over the epoch the caller
% actually supplied is the least surprising default; a narrower window stays
% available, but opt-in.
if ~any(strcmpi(fwd(1:2:end), 'detect_win_ms'))
    fwd(end+1:end+2) = {'detect_win_ms', [time_ms(1) time_ms(end)]};
end

srcMask = [];
if p.InteriorOnly, srcMask = interiorMask(gridMask); end

wantPhase = nargout >= 3;

needNull = ischar(p.RhoThr) || isstring(p.RhoThr);
SD = waveSourcePoints(cube, gridMask, time_ms, ...
        'SourceMask', srcMask, 'ReturnPhase', wantPhase, ...
        'NullDraws', ternary(needNull, p.NullDraws, 0), ...
        'rhoThr', ternary(needNull, 0.5, p.RhoThr), ...   % placeholder when auto
        'Verbose', false, fwd{:});

nul = SD.nullRho(isfinite(SD.nullRho));
if needNull
    if isempty(nul)
        error('waveEvents:noNull','Null calibration requested but no null samples returned.');
    end
    rhoThr = prctile(nul, p.NullPct);
else
    rhoThr = double(p.RhoThr);
end

% ---- assemble events (re-threshold at the calibrated value) ----
rows = {};
nEval = 0;
srcOf = zeros(0,2);          % (epoch, eval-point) of each event row, for the phase stack
for s = 1:numel(SD.rho)
    r = SD.rho{s}; if isempty(r), continue; end
    nEval = nEval + sum(isfinite(r));
    dd = [];
    if isfield(SD,'dirDeg') && numel(SD.dirDeg) >= s, dd = SD.dirDeg{s}; end
    hit = find(r > rhoThr);
    for j = hit(:).'
        dj = NaN; if numel(dd) >= j, dj = dd(j); end
        rows(end+1,:) = {s, time_ms(SD.ep_idx{s}(j)), r(j), ...
                         SD.source_rc{s}(j,1), SD.source_rc{s}(j,2), dj}; %#ok<AGROW>
        srcOf(end+1,:) = [s j]; %#ok<AGROW>
    end
end
vn = {'epoch','time_ms','rho','srcRow','srcCol','dirDeg'};
if isempty(rows)
    E = table('Size',[0 numel(vn)], 'VariableTypes',repmat({'double'},1,numel(vn)), ...
        'VariableNames',vn);
else
    E = cell2table(rows, 'VariableNames',vn);
end

% ---- phase pattern per event row (same order as E) ----
P = [];
if wantPhase
    P = nan(size(cube,2), size(cube,3), size(srcOf,1));
    for k = 1:size(srcOf,1)
        P(:,:,k) = SD.phase{srcOf(k,1)}(:,:,srcOf(k,2));
    end
end

nE = size(cube,1);
dw = SD.Params.detect_win_ms;
winSec = max(diff(dw),1)/1000;
S = struct('rhoThr',rhoThr, 'nullPct',p.NullPct, 'nullRho',nul, ...
           'nEpochs',nE, 'nEval',nEval, 'nEvents',height(E), ...
           'detectRate', height(E)/max(nEval,1), ...
           'eventsPerSec', height(E)/max(nE*winSec,eps), ...
           'epochsWithEvent', numel(unique(E.epoch))/max(nE,1));

if p.Verbose
    fprintf(['waveEvents: thr %.3f (%s) | %d events / %d eval pts (%.0f%%) | ' ...
             '%.1f events/s | %.0f%% of epochs\n'], ...
        rhoThr, ternary(needNull, sprintf('null p%g', p.NullPct), 'fixed'), ...
        S.nEvents, S.nEval, 100*S.detectRate, S.eventsPerSec, 100*S.epochsWithEvent);
    if ~isempty(nul)
        fprintf('   null rho: p50 %.3f p95 %.3f p99 %.3f (n=%d draws)\n', ...
            prctile(nul,50), prctile(nul,95), prctile(nul,99), numel(nul));
    end
end
end

% ---------------------------------------------------------------------
function m = interiorMask(gm)
gm = logical(gm);
up=[false(1,size(gm,2)); gm(1:end-1,:)]; dn=[gm(2:end,:); false(1,size(gm,2))];
lf=[false(size(gm,1),1) gm(:,1:end-1)];  rt=[gm(:,2:end) false(size(gm,1),1)];
m = gm & up & dn & lf & rt;
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end
