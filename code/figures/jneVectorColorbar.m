function cbax = jneVectorColorbar(fig, pos, cmap, cl, nLevels, labelStr, st, varargin)
% jneVectorColorbar(fig, pos, cmap, cl, nLevels, labelStr, st, ...)
% Draw a discrete colorbar as VECTOR filled rectangles (patches) in a new axes at
% normalized figure position pos = [x y w h]. Unlike MATLAB's colorbar (which
% exports as a pattern-wrapped raster image that Adobe Illustrator drops), this is
% pure vector and round-trips into Illustrator. Band colours match a contour plot's
% discrete colormap: cmap sampled to nLevels over [cl(1) cl(2)].
%
% OPTIONS: 'Ticks' (default [] = auto/round), 'TickLabels' (cellstr, matches Ticks).
q = inputParser;
q.addParameter('Ticks', []);
q.addParameter('TickLabels', {});
q.parse(varargin{:}); o = q.Results;

cbax = axes('Parent', fig, 'Position', pos);
edges = linspace(cl(1), cl(2), nLevels + 1);
cols  = cmap(round(linspace(1, size(cmap,1), nLevels)), :);
hold(cbax, 'on');
for k = 1:nLevels
    patch(cbax, [0 1 1 0], [edges(k) edges(k) edges(k+1) edges(k+1)], ...
        cols(k,:), 'EdgeColor', 'none');
end
hold(cbax, 'off');
set(cbax, 'XLim',[0 1], 'YLim',[cl(1) cl(2)], 'XTick',[], ...
    'YAxisLocation','right', 'Box','on', 'Layer','top', 'TickDir','out', ...
    'FontName',st.FontName, 'FontSize',st.FontSizeTick, 'LineWidth',st.LineWidthAxis, ...
    'TickLabelInterpreter','none');
if ~isempty(o.Ticks), set(cbax, 'YTick', o.Ticks); end
if ~isempty(o.TickLabels), set(cbax, 'YTickLabel', o.TickLabels); end
yl = ylabel(cbax, labelStr, 'FontName',st.FontName, 'FontSize',st.FontSizeLabel, ...
    'Interpreter','none');
yl.Rotation = 90;
end
