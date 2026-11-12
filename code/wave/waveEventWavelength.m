function [E, D] = waveEventWavelength(E, cube, gridMask, time_ms, varargin)
% waveEventWavelength   Per-event spatial wavelength, frequency and speed.
%
% PURPOSE
% -------
% waveEvents says WHEN a wave was detected and WHICH WAY it travelled; it says
% nothing about how long the wave was in space. This adds that, one row at a time,
% measured on the SAME phase field and at the SAME frame the detection was made on,
% so a wavelength always belongs to the event it is attached to.
%
% The phase field is rebuilt through waveAnalyticFrames -- the shared converter that
% reproduces the detector's 5-40 Hz + generalized_phase(gp_lp 0) construction. It is
% not re-implemented here, and by default a sample of events has its rho RECOMPUTED
% and checked against the stored value ('Verify'), so a silent divergence between
% this measurement and the detector becomes an error rather than a number.
%
% THE UNITS TRAP THAT MAKES THIS WORTH READING
% --------------------------------------------
% phase_gradient_complex_multiplication returns pm = |grad phi| / (2*pi). That is
% CYCLES per unit distance, not radians. So the wavelength is
%
%       lambda = 1 / k        with k = pm / pitch_mm   [cycles/mm]
%
% and NOT 2*pi/k. waveGP's meanGradMag is this same cycles-per-mm quantity; its
% header called it "rad/mm" until 2026-09-15. No published result moved, because
% every consumer used it only as a relative ranking, but anyone who wrote 2*pi/k
% would have been off by 6.28x. waveGP's own speed (|IF|/pm * pitch) is right, and
% is right precisely BECAUSE pm is in cycles -- the two agree here: speed = f*lambda.
%
% WHY THE AMPLITUDE WEIGHTING INVERTS ONLY ONCE, AT THE END
% --------------------------------------------------------
% k is averaged across electrodes and the result inverted, rather than averaging the
% per-electrode 1/k. Wherever the local gradient approaches zero -- which it does, on
% an array this much smaller than a wavelength -- 1/k diverges, and a mean of
% divergent terms is set by its worst electrode. Averaging in k-space is bounded and
% is the same order waveGP uses. The weight is the analytic amplitude with the
% bottom 'AmpQuantile' zeroed and fabricated (hole-filled) cells zeroed, identical
% to waveGP, so no interpolated cell enters a reported value.
%
% APERTURE -- READ BEFORE QUOTING AN ABSOLUTE NUMBER.
% On this array a typical lambda is several times the array's own diagonal, so it is
% an EXTRAPOLATION from a fraction of a cycle of phase, not a measurement of a cycle.
% It is heavy-tailed; measurement noise adds spurious gradient and so SHORTENS
% lambda, which means anything that improves SNR lengthens it mechanically; and its
% absolute value is not comparable across arrays of different extent.
% .lamOverAperture is returned for exactly this reason: use RANK statistics within an
% array, and treat any absolute lambda as an order of magnitude.
%
% INPUTS
% ------
% E         event table from waveEvents (needs .epoch, .time_ms; .rho .srcRow
%           .srcCol additionally required when 'Verify' is on)
% cube      nEpoch x nRow x nCol x T, the SAME cube the events came from
% gridMask  nRow x nCol logical
% time_ms   1 x T, the same axis passed to waveEvents
% Name/value:
%   'PixelSpacingMm' 0.7   electrode pitch; pass C.pitch_mm, do not accept this default
%   'Fs'             1000
%   'AmpQuantile'    0.20  matches waveGP's amplitude mask
%   'Verify'         25    recompute rho on this many events and assert a match
%                          (0 = off). Cheap, and it is the only thing tying this
%                          measurement to the detector.
%   'Verbose'        true
%
% OUTPUTS
% -------
% E  the input table with columns appended:
%      .k_cyc_mm        amplitude-weighted spatial frequency (cycles/mm)
%      .lambda_mm       1 / k_cyc_mm
%      .freqHz          amplitude-weighted |instantaneous frequency| (Hz)
%      .speed_mm_s      freqHz / k_cyc_mm  (== freqHz * lambda_mm)
%      .lamOverAperture lambda_mm / array diagonal -- greater than 1 means unresolved
%      .amp             array-mean analytic amplitude at the event frame (for amplitude-
%                       matching speeds across surfaces; added 2026-10-28)
% D  diagnostics: .apertureMm .pitchMm .nVerified .maxRhoErr .fracUnresolved
%
% EXAMPLE
%   C = saccadeEpochCube(S(k).folders);
%   [E,~] = waveEvents(C.cube, C.gridMask, C.time_ms, 'pixelSpacing_mm',C.pitch_mm);
%   E = waveEventWavelength(E, C.cube, C.gridMask, C.time_ms, 'PixelSpacingMm',C.pitch_mm);
%
% See also: waveEvents, waveSourcePoints, waveAnalyticFrames, waveGP

p = struct('PixelSpacingMm',0.7, 'Fs',1000, 'AmpQuantile',0.20, ...
           'Verify',25, 'Verbose',true);
for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end

gm = logical(gridMask);
pitch = p.PixelSpacingMm;
checkPitchDefault('waveEventWavelength', varargin, gm);

% Array aperture = diagonal of the bounding box of REAL electrodes. This is the
% length scale beyond which a wavelength is no longer being measured.
[rr, cc] = find(gm);
apertureMm = hypot((max(rr)-min(rr))*pitch, (max(cc)-min(cc))*pitch);

nEv = height(E);
k_cyc_mm = nan(nEv,1); freqHz = nan(nEv,1); ampEv = nan(nEv,1);
if nEv == 0
    E = addCols(E, k_cyc_mm, freqHz, apertureMm, ampEv);
    D = struct('apertureMm',apertureMm, 'pitchMm',pitch, 'nVerified',0, ...
               'maxRhoErr',NaN, 'fracUnresolved',NaN);
    return;
end

doVerify = p.Verify > 0 && ...
    all(ismember({'rho','srcRow','srcCol'}, E.Properties.VariableNames));
vSet = [];
if doVerify
    % PRIVATE stream: which events get checked must be reproducible WITHOUT
    % reseeding the caller's global RNG. An engine that calls rng() is how
    % waveTimecoursePermTest ended up resampling one subset repeatedly.
    rs = RandStream('twister','Seed',11);
    vSet = randperm(rs, nEv, min(p.Verify, nEv));
end
maxRhoErr = 0; nVerified = 0;

% One analytic-signal build per EPOCH, not per event: epochs carry many events and
% the GP construction is the expensive part.
epochs = unique(E.epoch(:))';
for e = epochs
    rows = find(E.epoch == e);
    fr = squeeze(cube(e,:,:,:));
    X  = waveAnalyticFrames(fr, gm);           % SHARED: fills holes, 5-40 Hz, GP

    pm = phase_gradient_complex_multiplication(X, 1);        % CYCLES per electrode
    K  = pm ./ pitch;                                        % cycles/mm
    IF = instantaneous_frequency(X, p.Fs);                   % Hz

    % amplitude weight, exactly as waveGP: sub-quantile zeroed, fabricated cells zeroed
    A = abs(X); T = size(A,3);
    thr = quantile(A(repmat(gm,1,1,T)), p.AmpQuantile);
    W = A; W(A < thr) = 0; W = W .* gm;

    for r = rows(:)'
        [~, kEv] = min(abs(time_ms - E.time_ms(r)));
        w  = W(:,:,kEv); sw = sum(w(:), 'omitnan');
        if ~(sw > 0), continue; end
        kf = K(:,:,kEv);  ff = abs(IF(:,:,kEv));
        k_cyc_mm(r) = sum(kf(:).*w(:), 'omitnan') / sw;
        freqHz(r)   = sum(ff(:).*w(:), 'omitnan') / sw;
        a = A(:,:,kEv); ampEv(r) = mean(a(gm), 'omitnan');    % for amplitude-matching speeds

        % Tie this measurement to the detector: same frame, same field, same rho.
        if ismember(r, vSet)
            phM = angle(X(:,:,kEv)); phM(~gm) = NaN;
            rhoChk = phase_correlation_distance(phM, [E.srcCol(r); E.srcRow(r)], pitch);
            maxRhoErr = max(maxRhoErr, abs(rhoChk - E.rho(r)));
            nVerified = nVerified + 1;
        end
    end
end

if nVerified > 0 && maxRhoErr > 1e-6
    error('waveEventWavelength:diverged', ...
        ['Recomputed rho differs from the stored event rho by %.3g on %d checked ' ...
         'events. The wavelength would not belong to the detection -- the cube, the ' ...
         'time axis or the phase construction does not match the one waveEvents used.'], ...
        maxRhoErr, nVerified);
end

E = addCols(E, k_cyc_mm, freqHz, apertureMm, ampEv);
D = struct('apertureMm',apertureMm, 'pitchMm',pitch, 'nVerified',nVerified, ...
           'maxRhoErr',maxRhoErr, ...
           'fracUnresolved',mean(E.lamOverAperture > 1, 'omitnan'));

if p.Verbose
    fprintf(['waveEventWavelength: %d events | aperture %.1f mm | lambda median %.1f mm ' ...
             '(%.1fx aperture) | %.0f%% unresolved | f median %.1f Hz | ' ...
             'verified %d (max rho err %.1e)\n'], ...
        nEv, apertureMm, median(E.lambda_mm,'omitnan'), ...
        median(E.lamOverAperture,'omitnan'), 100*D.fracUnresolved, ...
        median(E.freqHz,'omitnan'), nVerified, maxRhoErr);
end
end

% ======================================================================
function E = addCols(E, k_cyc_mm, freqHz, apertureMm, ampEv)
E.k_cyc_mm        = k_cyc_mm;
E.lambda_mm       = 1 ./ k_cyc_mm;
E.freqHz          = freqHz;
E.speed_mm_s      = freqHz ./ k_cyc_mm;
E.lamOverAperture = E.lambda_mm ./ apertureMm;
% array-mean analytic amplitude at the event frame. Phase-gradient read-outs track amplitude
% (k differs ~4x across amplitude quintiles; utilities/wave/README.md item 5), so any speed
% comparison between surfaces with different evoked amplitude must be matched on this.
E.amp             = ampEv;
end
