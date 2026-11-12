function cm = jneDisplayMap(cm)
% jneDisplayMap  The array orientation the JNE figures are DRAWN in.
%
% The figures use the channel map the Kuzum group supplied with the array (row 1 =
% channels 57 63 1 7). The lab pipeline (getChannelMap) uses the anatomical orientation
% adopted 2026-08-05, which is that map turned upside down (row 1 = 53 33 31 11). So the
% figures flip rows for display only: every number is computed on the pipeline map and is
% unchanged by this (Blake, 2026-11-12). Used by jneData (Fig 7) and jneFlipEpoch (Fig 8a,
% Supplementary Videos).
cm = flipud(cm);
end
