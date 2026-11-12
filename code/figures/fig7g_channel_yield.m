function fig = fig7g_channel_yield(D, st, varargin)
% FIG 7g (JNE) — Functional channel yield (% responsive contacts) vs implant age.
% Shows five representative timepoints spanning the recording period
% (1 wk / 2 wk / ~7 wk / 2 mo / 4 mo). Uniform gray line + markers (no per-session
% highlighting). Set 'Dates',[] to instead plot all sessions <= 'MaxDays'
% (excluding 'ExcludeDates').
if nargin < 1 || isempty(D),  D  = jneData();  end
if nargin < 2 || isempty(st), st = jneStyle(); end
p = inputParser;
p.addParameter('YLim', [50 100], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Dates', ["251024" "251028" "251201" "251219" "260203"], @(x)isstring(x)||ischar(x)||iscellstr(x)); % 1wk/2wk/~7wk(d48)/2mo/4mo
p.addParameter('MaxDays', 130, @isscalar);              % used only when Dates is empty
p.addParameter('ExcludeDates', "251211", @(x)isstring(x)||ischar(x)||iscellstr(x));
p.addParameter('Style', 'line', @(s)ischar(s)||isstring(s));   % 'line' | 'stem'
p.addParameter('SaveTag', 'fig7g_channel_yield', @(s)ischar(s)||isstring(s));
p.addParameter('SaveDir', jneOutDir(), @(s)ischar(s)||isstring(s));
p.addParameter('Save', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt = p.Results;

Q = D.Q;
if isempty(opt.Dates)
    keep = Q.ageDays <= opt.MaxDays & ~ismember(string(Q.Date), string(opt.ExcludeDates));
else
    keep = ismember(string(Q.Date), string(opt.Dates));
end
Q = Q(keep,:);
[~, o] = sort(Q.ageDays); Q = Q(o,:);         % chronological
x = Q.ageDays; y = Q.respFrac*100;

fig = jneFigure(st.WidthSingle, st.RowHeight, st);
ax = axes(fig); hold(ax,'on');
if strcmpi(opt.Style,'stem')
    stem(ax, x, y, 'filled', 'Color', st.Gray, 'MarkerFaceColor', st.Black, ...
        'MarkerEdgeColor', st.Black, 'MarkerSize', 4, 'LineWidth', st.LineWidthData, ...
        'BaseValue', opt.YLim(1));
else
    plot(ax, x, y, '-', 'Color', st.Gray, 'LineWidth', st.LineWidthAxis);
    scatter(ax, x, y, 24, st.Black, 'filled');      % uniform black dots, no rep highlight
end
hold(ax,'off');
ylim(ax, opt.YLim); xlim(ax, [0 max(x)+10]);
xlabel(ax, 'Days post-implant'); ylabel(ax, 'Functional channels (%)');
title(ax, 'Stimulus-responsive channels over time');
jneApplyAxes(ax, st);

if opt.Save, jneSave(fig, opt.SaveDir, char(opt.SaveTag), st); end
end
