function jneWaveGrid(ax, erp, tMs, channelMap, ampScale, col, xlMs, st)
% jneWaveGrid(ax, erp, tMs, channelMap, ampScale, col, xlMs, st)
% Draw a small trial-averaged ERP trace at each contact's array (row,col)
% location over the display window xlMs. ampScale (µV) sets the shared vertical
% scaling so multiple panels are directly comparable.
[nr, nc] = size(channelMap);
chMap1 = channelMap + 1;
kt = tMs >= xlMs(1) & tMs <= xlMs(2);
tw = tMs(kt);
tn = (tw - tw(1)) / (tw(end) - tw(1));
cellW = 0.82; cellH = 0.42;
hold(ax,'on');
for rr = 1:nr
    for cc = 1:nc
        ch = chMap1(rr,cc);
        if isnan(ch) || ch<1 || ch>size(erp,1), continue; end
        x = cc + (tn - 0.5)*cellW;
        y = (nr - rr + 1) + (erp(ch,kt)/ampScale)*cellH;
        plot(ax, x, y, 'Color', col, 'LineWidth', st.LineWidthThin);
    end
end
% scale bar, tucked in the empty band BELOW the array (clear of all traces):
% L-shaped, amplitude label to the LEFT of the vertical stroke (horizontal, right
% aligned), time label BELOW the horizontal stroke — well separated, no overlap.
xb = 0.9; yb = -0.45;
xs = (100/(tw(end)-tw(1)))*cellW;
plot(ax, [xb xb], [yb yb+cellH], 'k-', 'LineWidth', st.LineWidthAxis);   % amplitude
plot(ax, [xb xb+xs], [yb yb], 'k-', 'LineWidth', st.LineWidthAxis);      % 100 ms
text(ax, xb-0.16, yb+cellH/2, sprintf('%.0f %s', ampScale, st.muV), ...
    'HorizontalAlignment','right','VerticalAlignment','middle','Clipping','off', ...
    'Interpreter','none','FontSize',st.FontSizeAnnot,'FontName',st.FontName);
text(ax, xb+xs/2, yb-0.16, '100 ms', 'HorizontalAlignment','center', ...
    'VerticalAlignment','top','Clipping','off', ...
    'Interpreter','none','FontSize',st.FontSizeAnnot,'FontName',st.FontName);
hold(ax,'off');
axis(ax,'equal'); axis(ax,'off');
xlim(ax, [-1.1 nc+0.7]); ylim(ax, [-1.0 nr+1]);
end
