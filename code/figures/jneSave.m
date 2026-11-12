function jneSave(fig, outDir, tag, st)
% jneSave(fig, outDir, tag, st) — export a figure uniformly for the JNE set:
% a raster PNG (preview) and a vectorized SVG (Illustrator-editable text).
if nargin < 4 || isempty(st), st = jneStyle(); end
if ~isfolder(outDir), mkdir(outDir); end
base = fullfile(outDir, char(tag));
t0 = now;
try
    exportgraphics(fig, [base '.png'], 'Resolution', st.ExportDPI);
catch ME
    warning('jneSave:png', 'PNG export failed: %s', ME.message);
end
% A PNG that is open elsewhere (an image viewer, a preview pane) can be LOCKED on Windows;
% the export then fails and the OLD file stays on disk looking current (2026-10-28). Say so.
dp = dir([base '.png']);
if isempty(dp) || dp.datenum < t0 - 1/86400
    fprintf(2, '[jneSave] *** %s.png was NOT updated (file locked or export failed) -- it is STALE ***\n', tag);
end
% SVG via exportgraphics vector: carries colorbars (which print -dsvg can drop on
% imagesc panels) and, with Interpreter='none' text, stays editable in Illustrator.
try
    exportgraphics(fig, [base '.svg'], 'ContentType', 'vector');
    fixSvgBase64([base '.svg']);   % strip whitespace after 'base64,' (Illustrator won't decode it)
catch ME
    warning('jneSave:svg', 'SVG export failed (%s); falling back to print', ME.message);
    try, print(fig, [base '.svg'], '-dsvg', '-vector'); fixSvgBase64([base '.svg']); catch, end
end
end


function fixSvgBase64(f)
% exportgraphics writes embedded rasters as "data:image/png;base64, <data>" — the
% space after the comma makes Adobe Illustrator silently drop the image (heatmaps
% AND colorbars vanish). Remove it. Byte-level so UTF-8 glyphs (µ, ²) are preserved.
try
    fid = fopen(f, 'r'); if fid < 0, return; end
    b = fread(fid, Inf, '*uint8')'; fclose(fid);
    s = char(b);                                 % 1 byte <-> 1 char code (no UTF-8 decode)
    s2 = strrep(s, ';base64, ', ';base64,');
    if ~strcmp(s, s2)
        fid = fopen(f, 'w'); if fid < 0, return; end
        fwrite(fid, uint8(s2), 'uint8'); fclose(fid);
    end
catch
end
end
