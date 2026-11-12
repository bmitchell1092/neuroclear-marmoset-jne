function [p2p, latMs] = jneWindowMetrics(erp, tMs, win)
% [p2p, latMs] = jneWindowMetrics(erp, tMs, win)
% Per-channel peak-to-peak amplitude and peak latency of the trial-averaged ERP
% within the window win = [lo hi] ms. erp is nChan x nT.
sel = tMs >= win(1) & tMs <= win(2);
E  = erp(:, sel);
tt = tMs(sel);
p2p = max(E, [], 2) - min(E, [], 2);
[~, imx] = max(abs(E), [], 2);       % time of largest |deflection|
latMs = tt(imx);
latMs = latMs(:);
end
