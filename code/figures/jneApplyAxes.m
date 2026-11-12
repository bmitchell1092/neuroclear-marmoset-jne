function jneApplyAxes(ax, st)
% jneApplyAxes(ax, st) — stamp the JNE style onto an existing axes (fonts, line
% width, box/ticks) and bump its label/title font sizes. Call after building the
% axes content and setting xlabel/ylabel/title.
if nargin < 2 || isempty(st), st = jneStyle(); end
set(ax, 'FontName',st.FontName, 'FontSize',st.FontSizeTick, ...
    'LineWidth',st.LineWidthAxis, 'Box',st.Box, 'TickDir',st.TickDir, ...
    'TickLength',st.TickLength, 'Layer','top', 'GridAlpha',st.GridAlpha, ...
    'TickLabelInterpreter','none');   % 'none' -> SVG text stays editable in Illustrator
if ~isempty(ax.XLabel), ax.XLabel.FontName = st.FontName; ax.XLabel.FontSize = st.FontSizeLabel; ax.XLabel.Interpreter = 'none'; end
if ~isempty(ax.YLabel), ax.YLabel.FontName = st.FontName; ax.YLabel.FontSize = st.FontSizeLabel; ax.YLabel.Interpreter = 'none'; end
if ~isempty(ax.Title),  ax.Title.FontName  = st.FontName; ax.Title.FontSize  = st.FontSizeTitle; ax.Title.Interpreter = 'none'; end
end
