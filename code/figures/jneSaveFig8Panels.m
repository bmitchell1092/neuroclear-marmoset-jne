function jneSaveFig8Panels(fig, saveDir, tag, st, H, Hd)
% jneSaveFig8Panels  Save each panel of Fig 8 as its own file (<tag>_a, _b, _c[, _d, _e]).
%
% Everything outside a panel is hidden, the figure is exported (exportgraphics crops to what
% is visible), and visibility is restored -- so the panel files are exactly what is in the
% figure (same sizes, fonts, data), not a redraw. Panels are told apart by position, using
% the fig8_travelingWaves layout: a = the frame grid on the left (incl. its colourbars and
% headers), b = top right, c = bottom right, and with the longevity row, d/e below it.
% H, Hd: heights (inches) of the a-c block and of the longevity row (0 when off).
Ht = H + Hd; yRow = Hd / Ht;                         % top of the longevity row
yMid = (0.5*H + Hd) / Ht;                            % between b and c
letters = 'abc'; if Hd > 0, letters = 'abcde'; end
h = [findall(fig, 'Type','axes'); findall(fig, 'Type','textboxshape')];
P = cell2mat(get(h, 'Position'));
low = Hd > 0 & P(:,2) < yRow;
sel = struct('a', ~low & P(:,1) < 0.67, 'b', ~low & P(:,1) >= 0.67 & P(:,2) >= yMid, ...
             'c', ~low & P(:,1) >= 0.67 & P(:,2) < yMid, 'd', low & P(:,1) < 0.48, ...
             'e', low & P(:,1) >= 0.48);
% every object to toggle: each axes with all its contents (lines, patches, titles, labels),
% and each annotation; remember the original state
parts = cell(numel(h), 1);
for i = 1:numel(h), parts{i} = findall(h(i)); end
all0 = vertcat(parts{:}); vis0 = get(all0, 'Visible');
for L = letters
    for i = 1:numel(h)
        if sel.(L)(i), set(parts{i}, {'Visible'}, vis0(ismember(all0, parts{i})));
        else, set(parts{i}, 'Visible', 'off'); end
    end
    jneSave(fig, saveDir, sprintf('%s_%s', tag, L), st);
end
set(all0, {'Visible'}, vis0);                         % restore
end
