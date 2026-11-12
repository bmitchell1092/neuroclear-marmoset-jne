function fig = fig7f_psd_progression(D, st, varargin)
% FIG 7f (JNE) — Resting (spontaneous) power spectrum at 2 wk / 2 mo / 4 mo, one
% median trace per session (colored by age), showing spectral stability.
if nargin < 1 || isempty(D),  D  = jneData();  end
if nargin < 2 || isempty(st), st = jneStyle(); end
p = inputParser;
p.addParameter('FreqLim', [1 100], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('NotchHz', [30 60 90 120 180], @isnumeric);   % line noise + harmonics
p.addParameter('NotchBW', 2, @isscalar);                     % half-width (Hz) per notch
p.addParameter('GrayRange', [0 0.6], @(x)isnumeric(x)&&numel(x)==2); % black->gray gradient by age
p.addParameter('Alpha', 0.8, @isscalar);                     % trace transparency (show overlaps)
p.addParameter('SaveDir', jneOutDir(), @(s)ischar(s)||isstring(s));
p.addParameter('Save', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt = p.Results;

nP = numel(D.darkPSD);
grays = linspace(opt.GrayRange(1), opt.GrayRange(2), max(nP,2));   % 2 wk black -> 4 mo gray
fig = jneFigure(st.WidthSingle, st.RowHeight, st);
ax = axes(fig); hold(ax,'on');
h = gobjects(0); labs = strings(0);
for i = 1:nP
    dp = D.darkPSD(i);
    if isempty(dp.f), continue; end
    med = median(dp.Pxx, 1, 'omitnan');
    med = notchInterp(dp.f, med, opt.NotchHz, opt.NotchBW);   % interp over line-noise bins
    hi = plot(ax, dp.f, 10*log10(med+eps), 'Color', grays(i)*[1 1 1], ...
        'LineWidth', st.LineWidthData);
    hi.Color(4) = opt.Alpha;                                 % semi-transparent so overlaps show
    h(end+1) = hi; labs(end+1) = dp.label; %#ok<AGROW>
end
hold(ax,'off');
set(ax, 'XScale', 'log'); xlim(ax, opt.FreqLim);
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, st.psdYLabel);
title(ax, 'Resting spectrum over time');
if ~isempty(h)
    lg = legend(ax, h, cellstr(labs), 'Location','southwest'); lg.Box = 'off';
    lg.FontName = st.FontName; lg.FontSize = st.FontSizeAnnot; lg.Interpreter = 'none';
end
jneApplyAxes(ax, st);

if opt.Save, jneSave(fig, opt.SaveDir, 'fig7f_psd_progression', st); end
end


function p = notchInterp(f, p, notchHz, bw)
% Remove narrow line-noise spikes for display, leaving the broadband spectrum
% (incl. the broad alpha peak) untouched. Two passes on log-power:
%   1) explicit notch bands (line frequency + harmonics), +/- bw Hz.
%   2) automatic narrow-spike detection: bins that exceed a local moving-median
%      baseline by > thr (log10 units ~ 3.5 dB). movmedian tracks broad features
%      (alpha) so only sharp artifacts are flagged.
% Flagged bins are then filled by linear interpolation from the surviving bins.
f = f(:).'; p = p(:).'; lp = log10(p + eps); n = numel(f);
bad = false(1, n);
for fc = notchHz(:).'
    bad = bad | (f >= (fc-bw) & f <= (fc+bw));
end
win = max(5, 2*round(2/median(diff(f)))+1);      % ~4 Hz moving-median window (odd)
base = movmedian(lp, win, 'omitnan');
bad = bad | ((lp - base) > 0.35);
bad(1) = false; bad(end) = false;                % keep endpoints as interp anchors
good = ~bad;
if any(bad) && nnz(good) >= 2
    lp(bad) = interp1(f(good), lp(good), f(bad), 'linear');
end
p = 10.^lp;
end
