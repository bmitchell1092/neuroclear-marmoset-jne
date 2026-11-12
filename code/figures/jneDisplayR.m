function R = jneDisplayR(R)
% jneDisplayR  A flashERP result (or its slim jneData copy) in the figure orientation:
% channelMap and every R.map grid flipped by jneDisplayMap. Per-channel vectors (erp,
% ampP2P, responsive, ...) are indexed by channel and need no change. Display only.
R.channelMap = jneDisplayMap(R.channelMap);
if isfield(R, 'map') && isstruct(R.map)
    for f = fieldnames(R.map)'
        if ~isempty(R.map.(f{1})), R.map.(f{1}) = jneDisplayMap(R.map.(f{1})); end
    end
end
end
