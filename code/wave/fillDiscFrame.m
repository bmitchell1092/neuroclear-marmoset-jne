function Fi = fillDiscFrame(F, method)
% fillDiscFrame   Fill non-electrode / bad-channel holes for wave analysis, safely.
%
% Single source of truth for the hole-filling used everywhere in the traveling-
% wave code (disc corners + interpolated bad channels). Interpolates the missing
% cells, then CLAMPS the fabricated values to the observed real-electrode range so
% a hole on the convex-hull boundary (e.g. a cluster of bad channels at the array
% edge) can never extrapolate to an out-of-range spike. A nearest-neighbour pass
% backfills anything the interpolant left as NaN, so no NaN survives for the phase
% gradient.
%
% Rationale: qcInterpolation found that plain fillmissing2('linear') blows up at
% edge bad-clusters (a top-row 2x2 of bad channels became ~-15 uV, far outside its
% green neighbours). Clamping to the observed range keeps smooth interior gradients
% while bounding the boundary, giving 0% overshoot without the blockiness of pure
% 'nearest'. See utilities/wave/qcInterpolation.m.
%
% INPUT  F       nRow x nCol frame with NaN at cells to fill.
%        method  fillmissing2 method for the interior (default 'linear').
% OUTPUT Fi      fully filled frame, fabricated cells clamped to range(F).
%
% See also: fillmissing2, waveGP, waveSourcePoints, grandAverageWaveMovie
if nargin < 2, method = 'linear'; end
if ~any(isnan(F(:))), Fi = F; return; end
lo = min(F(:),[],'omitnan'); hi = max(F(:),[],'omitnan');
Fi = fillmissing2(F, method);
Fi = min(max(Fi, lo), hi);                       % clamp fabricated values to observed range
if any(isnan(Fi(:))), Fi = fillmissing2(Fi,'nearest'); end   % guarantee no NaN for the gradient
end
