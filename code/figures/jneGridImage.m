function jneGridImage(ax, vals, channelMap, cmap, cl, st, varargin)
% jneGridImage(ax, vals, channelMap, cmap, cl, st, ...)
% Draw a per-channel value vector laid on the electrode grid (channelMap,
% 0-indexed) as either a smooth interpolated image or a few-level filled contour.
% Missing/absent contacts are transparent (show the axes NaN color).
%
% OPTIONS (name/value)
%   'UpN'          fine-grid resolution (default 220)
%   'Mode'         'smooth' (imagesc) | 'contour' (filled contour bands) (default 'smooth')
%   'NLevels'      number of contour bands for 'contour' (default 5)
%   'Interp'       'natural' (smoothest) | 'linear' (tighter/more local) (default 'natural')
%   'SmoothSigma'  extra Gaussian blur on the interpolated grid, fine-grid px (default 0)
%   'ShowContacts' overlay a marker at each real contact center (default false)
%   'ContactAlpha' contact marker face alpha (default 0.30)
%   'ContactSize'  contact marker area (default 26)
q = inputParser;
q.addParameter('UpN', 220);
q.addParameter('Mode', 'smooth');
q.addParameter('NLevels', 5);
q.addParameter('Interp', 'natural');
q.addParameter('SmoothSigma', 0);
q.addParameter('ShowContacts', false);
q.addParameter('ContactAlpha', 0.30);
q.addParameter('ContactSize', 26);
q.parse(varargin{:}); o = q.Results;

[nr, nc] = size(channelMap);
chMap1 = channelMap + 1;
xs = []; ys = []; vs = [];
for rr = 1:nr
    for cc = 1:nc
        ch = chMap1(rr,cc);
        if ~isnan(ch) && ch>=1 && ch<=numel(vals) && isfinite(vals(ch))
            xs(end+1)=cc; ys(end+1)=rr; vs(end+1)=vals(ch); %#ok<AGROW>
        end
    end
end

xf = linspace(1, nc, o.UpN);
yf = linspace(1, nr, o.UpN);
[XF, YF] = meshgrid(xf, yf);
if numel(vs) >= 3
    F = scatteredInterpolant(xs(:), ys(:), vs(:), o.Interp, 'none');
    VF = F(XF, YF);
else
    VF = nan(size(XF));
end
if o.SmoothSigma > 0, VF = gaussSmoothNaN(VF, o.SmoothSigma); end
if isempty(cl), cl = [min(vs) max(vs)]; end

hold(ax, 'on');
switch lower(string(o.Mode))
    case "contour"
        edges = linspace(cl(1), cl(2), o.NLevels+1);
        VFc = VF;                              % clamp finite values into [cl] so contourf fills
        VFc(VFc < cl(1)) = cl(1);              % beyond-range regions would otherwise be left
        VFc(VFc > cl(2)) = cl(2);              % blank (background shows through -> looks missing)
        contourf(ax, XF, YF, VFc, edges, 'LineColor',[0.35 0.35 0.35], 'LineWidth',0.4);
        colormap(ax, cmap(round(linspace(1, size(cmap,1), o.NLevels)), :));  % N discrete bands
    otherwise
        % Smooth raster map (fine for PNG/screen). NOTE: MATLAB exports imagesc as a
        % <pattern> fill that Adobe Illustrator can't round-trip — use 'contour' mode
        % for SVGs destined for Illustrator (vector patches, editable).
        h = imagesc(ax, xf, yf, VF);
        set(h, 'AlphaData', isfinite(VF));
        colormap(ax, cmap);
end
if o.ShowContacts
    scatter(ax, xs, ys, o.ContactSize, 'o', 'MarkerFaceColor',[1 1 1], ...
        'MarkerFaceAlpha',o.ContactAlpha, 'MarkerEdgeColor',[0.15 0.15 0.15], ...
        'MarkerEdgeAlpha',0.6, 'LineWidth',0.4);
end
hold(ax, 'off');

set(ax, 'YDir','reverse', 'Color', st.NaNColor);
axis(ax, 'image');
set(ax, 'CLim', cl, 'XTick', [], 'YTick', []);
xlim(ax, [0.5 nc+0.5]); ylim(ax, [0.5 nr+0.5]);
end


% =====================================================================
function Vs = gaussSmoothNaN(V, sigma)
% Gaussian smooth that ignores NaNs (normalized convolution), preserving the
% NaN (outside-hull) region.
m = isfinite(V);
V0 = V; V0(~m) = 0;
num = imgaussfilt(V0, sigma);
den = imgaussfilt(double(m), sigma);
Vs = num ./ den;
Vs(~m) = NaN;
end
