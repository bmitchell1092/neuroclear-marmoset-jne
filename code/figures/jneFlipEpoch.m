function P = jneFlipEpoch(P)
% jneFlipEpoch  Turn a jneSpontEpochData / jneFlashEpochData result into the figure
% orientation (jneDisplayMap: rows flipped). Display only -- detection ran on the pipeline
% grid and is not redone.
%
% Frames (Z, ph, gm, cube, cubePadded) are flipped top-to-bottom. Each detected wave is
% mirrored with them: its direction (top-CW degrees) a -> 180 - a, so its vertical
% component reverses, and its source row r -> nRows + 1 - r.
nr = size(P.gm, 1);
P.gm = jneDisplayMap(P.gm);
for f = ["Z","ph"]
    if isfield(P, f), P.(f) = flip(P.(f), 1); end
end
for f = ["cube","cubePadded"]                       % trial x row x col x time
    if isfield(P, f), P.(f) = flip(P.(f), 2); end
end
if isfield(P, 'E') && height(P.E) > 0
    P.E.dirDeg = mod(180 - P.E.dirDeg + 180, 360) - 180;
    P.E.srcRow = nr + 1 - P.E.srcRow;
end
P = rmfield(P, intersect(fieldnames(P), {'arrowDir','arrowCons'}));   % not drawn; not flipped
P.displayMap = 'jneDisplayMap (group-supplied orientation)';
end
