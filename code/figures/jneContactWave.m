function jneContactWave(ax, R, ch, xl, st)
% jneContactWave(ax, R, ch, xl, st)
% Draw ONE contact's trial-averaged flash VEP (mean +/- SEM) into ax, y-scaled to
% its own peak-to-peak so the waveform SHAPE fills the box regardless of absolute
% amplitude. Box on, no ticks, dashed stim-onset guide at t=0, peak-to-peak
% amplitude labelled top-left. Shared by the fig7d SEM insets and the standalone
% jne_contactInset figures so both look identical.
%
%   R   : slim flashERP result (needs .erp, .tMs, optional .sem)
%   ch  : erp-row index (1-based)
%   xl  : [lo hi] ms display window
%   st  : jneStyle struct
k  = R.tMs >= xl(1) & R.tMs <= xl(2);
t  = R.tMs(k);
mu = R.erp(ch, k);
if isfield(R,'sem') && ~isempty(R.sem), se = R.sem(ch, k); else, se = zeros(size(mu)); end
hold(ax,'on');
lo = min(mu-se); hi = max(mu+se); pad = 0.10*(hi-lo+eps);
if xl(1) < 0 && xl(2) > 0
    plot(ax, [0 0], [lo-pad hi+pad], ':', 'Color',[0.6 0.6 0.6], 'LineWidth',0.5);
end
patch(ax, [t fliplr(t)], [mu-se fliplr(mu+se)], st.Gray, 'FaceAlpha',0.35, 'EdgeColor','none');
plot(ax, t, mu, 'Color', st.Black, 'LineWidth', st.LineWidthData);
hold(ax,'off');
xlim(ax, xl); ylim(ax, [lo-pad hi+pad]);
set(ax, 'FontName',st.FontName, 'FontSize',st.FontSizeAnnot-1, ...
    'Box','on', 'LineWidth',st.LineWidthAxis, 'XTick',[], 'YTick',[]);
text(ax, 0.04, 0.96, sprintf('%.0f %s', max(mu)-min(mu), st.muV), 'Units','normalized', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', ...
    'FontName',st.FontName, 'FontSize',st.FontSizeAnnot-1);
end
