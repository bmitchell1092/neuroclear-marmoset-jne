function [dirDeg, cons] = waveMeanArrow(X, U, gm, k)
% waveMeanArrow   The ONE arrow drawn on a wave exemplar frame.
%
% PURPOSE
% -------
% Exemplar strips show a single centred arrow per frame, never a per-electrode vector field
% (readers found 64 arrows too busy to follow; Blake's convention since 2026-09-08). This is
% that arrow's definition, in one place, so every renderer that draws it -- the analysis
% montage (waveExemplarMontage 'phase+mean'), the movie, and the manuscript figures -- shows
% the same quantity.
%
%   direction    amplitude-weighted mean of the per-electrode propagation vectors
%   consistency  |mean of the UNIT vectors| in [0,1]: how much the array agrees; the
%                renderer scales arrow length by it, so a disorganised frame draws short
%
% INPUTS
%   X   nRow x nCol x T analytic signal (waveAnalyticFrames)
%   U   nRow x nCol x T propagation field (propagationField(X, gm, Fs))
%   gm  nRow x nCol logical grid mask
%   k   frame index (scalar or vector)
%
% OUTPUTS
%   dirDeg  screen frame (0 deg = right, CCW+) -- the frame of the array image, so it can be
%           drawn directly: dx = cosd(dirDeg), dy = -sind(dirDeg) on a YDir='reverse' axis
%   cons    directional consistency in [0,1]; NaN where the frame has no usable field
%
% See also: waveExemplarMontage, propagationField, waveAnalyticFrames

dirDeg = nan(size(k)); cons = nan(size(k));
for q = 1:numel(k)
    uu = U(:,:,k(q)); aa = abs(X(:,:,k(q)));
    uu(~gm) = NaN; aa(~gm) = NaN;
    gs = sum(uu(:).*aa(:), 'omitnan') / (sum(aa(:), 'omitnan') + eps);
    un = uu(isfinite(uu) & abs(uu) > 0);
    if ~isempty(un) && isfinite(gs) && abs(gs) > 0
        dirDeg(q) = rad2deg(angle(gs));
        cons(q)   = abs(mean(un ./ abs(un)));
    end
end
end
