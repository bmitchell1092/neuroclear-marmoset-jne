function X = waveAnalyticFrames(fr, gm)
% waveAnalyticFrames   Cube frames -> the analytic signal the DETECTOR sees.
%
% PURPOSE
% -------
% One place where "LFP frames" becomes "generalized phase", so every viewer, ranker and
% movie shows the same signal the statistics were computed on.
%
% ⚠ THIS IS THE BUG THAT KEEPS COMING BACK. `waveSourcePoints` bandpasses 5-40 Hz BEFORE
% generalized_phase (gp_lp = 0, matching the Muller demo). Feeding raw broadband LFP_filt
% straight to GP produces a far patchier field than the detector ever operates on: it
% makes real waves look like noise, and it made the movies appear to run ~6x faster than
% the analysis (measured 90.9 Hz implied vs the detector's 13.9 Hz). That mistake was
% found in FOUR separate renderers on 2026-08-31. Anything that displays or scores wave
% phase calls this function; nothing re-implements it.
%
% Holes are filled ONLY so the phase gradient is continuous. The caller re-masks with gm
% for display, so nothing fabricated is ever shown as an electrode.
%
% INPUTS
%   fr  nRow x nCol x T frames for ONE epoch
%   gm  nRow x nCol logical grid mask (unused in the maths; kept so callers pass intent)
%
% OUTPUT
%   X   nRow x nCol x T complex analytic signal (angle(X) = generalized phase)
%
% See also: waveExemplarMontage, waveRankExemplars, waveSourcePoints, waveEpochMovie

for k = 1:size(fr,3)
    f = fr(:,:,k);
    % fillDiscFrame(F, method) -- the 2nd arg is an interpolation METHOD, not a mask.
    if any(isnan(f(:))), fr(:,:,k) = fillDiscFrame(f); end
end

[nr,nc,T] = size(fr);
if exist('bandpass_filter','file') == 2
    fr = bandpass_filter(double(fr), 5, 40, 4, 1000);
else
    warning('waveAnalyticFrames:band', ...
        'bandpass_filter not on path; phase will look noisier than the detector sees');
end
if exist('generalized_phase','file') == 2
    X = generalized_phase(fr, 1000, 0);   % takes the nRow x nCol x T cube directly
else
    warning('waveAnalyticFrames:gp', ...
        'generalized_phase not on path; using hilbert (display only)');
    X = reshape(hilbert(reshape(fr,[nr*nc,T]).').', [nr,nc,T]);
end
end
