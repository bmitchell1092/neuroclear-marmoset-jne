function fig = fig7c_polarity_map(D, st, varargin)
% FIG 7c (JNE) — Representative early session (1 wk) spatial map of the flash-
% evoked potential. Default 'Metric' is 'polarity': the SIGNED evoked voltage at
% the population-peak latency on a blue-white-red diverging map, showing the
% spatial polarity (source/sink) / dipole structure of the evoked response.
% 'Metric','p2p' instead maps peak-to-peak transient amplitude (turbo).
%
% 'Date' selects the session (default '251024' = 1 wk). If that date is one of
% the cached representatives it is reused; otherwise its flashERP result is
% taken from D.other (jneData).
if nargin < 1 || isempty(D),  D  = jneData();  end
if nargin < 2 || isempty(st), st = jneStyle(); end
p = inputParser;
p.addParameter('Date', '251024', @(s)ischar(s)||isstring(s));      % 1 wk representative
p.addParameter('Metric', 'polarity', @(s)ischar(s)||isstring(s));  % 'polarity' | 'p2p'
p.addParameter('TransientMs', [50 100], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('PolarityWinMs', [40 150], @(x)isnumeric(x)&&numel(x)==2); % search window for the peak latency
p.addParameter('SnapshotMs', 130, @(x)isempty(x)||isscalar(x));    % fixed snapshot latency (ms); []=auto population-peak
p.addParameter('CLim', [-30 30], @(x)isempty(x)||numel(x)==2);     % colorbar range (µV); p2p coerced to non-negative
p.addParameter('Mode', 'contour', @(s)ischar(s)||isstring(s));     % 'smooth' | 'contour'
p.addParameter('NLevels', 10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('InterpBadChans', true, @(x)islogical(x)&&isscalar(x)); % fill QC-bad contacts from good neighbours
p.addParameter('BadChans', [], @isnumeric);                       % extra erp-rows to drop/interpolate over
p.addParameter('Interp', 'natural', @(s)ischar(s)||isstring(s));
p.addParameter('SmoothSigma', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ShowContacts', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('SaveTag', '', @(s)ischar(s)||isstring(s));
p.addParameter('SaveDir', jneOutDir(), @(s)ischar(s)||isstring(s));
p.addParameter('Save', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt = p.Results;

[R, lab, goodMask] = resolveSession(D, string(opt.Date));
sel = R.responsive; if ~any(sel), sel = true(R.nChan,1); end

if strcmpi(opt.Metric,'polarity')
    % signed voltage at a fixed snapshot latency, or the peak population deflection
    if ~isempty(opt.SnapshotMs)
        [~, it] = min(abs(R.tMs - opt.SnapshotMs));
    else
        kk = R.tMs >= opt.PolarityWinMs(1) & R.tMs <= opt.PolarityWinMs(2);
        ga = mean(abs(R.erp(sel,:)), 1, 'omitnan'); ga(~kk) = -inf;
        [~, it] = max(ga);
    end
    tStar = R.tMs(it);
    vals = R.erp(:, it);
    cmap = st.CMap_div;
    if isempty(opt.CLim) || opt.CLim(1) >= 0                 % need symmetric range
        M = prctile(abs(vals(sel)), 98); if ~isfinite(M)||M<=0, M = max(abs(vals(:))); end
        cl = [-M M];
    else
        cl = opt.CLim;
    end
    cbLabel = ['Evoked amplitude (' st.muV ')'];
    titleStr = sprintf('Evoked polarity (%d ms), %s', round(tStar), lab);
    defTag = 'fig7c_polarity_map';
else
    vals = jneWindowMetrics(R.erp, R.tMs, opt.TransientMs);
    cmap = st.CMap_p2p;
    cl = opt.CLim; if isempty(cl), cl = []; elseif cl(1) < 0, cl = [0 max(abs(cl))]; end  % p2p non-negative
    cbLabel = ['Peak-to-peak (' st.muV ')'];
    titleStr = sprintf('Evoked p2p (%d-%d ms), %s', opt.TransientMs(1), opt.TransientMs(2), lab);
    defTag = 'fig7c_p2p_map';
end

% Optionally interpolate over QC-bad / listed contacts. Their values are FILLED
% from good neighbours (grid positions kept) rather than dropped, so the array
% footprint/shape is preserved (dropping edge contacts would shrink the hull).
drop = false(numel(vals),1);
if opt.InterpBadChans, drop = drop | ~goodMask(:); end
if ~isempty(opt.BadChans), drop(opt.BadChans) = true; end
if any(drop)
    vals = jneFillContacts(vals, R.channelMap, drop);
    fprintf('[fig7c] filled %d contact(s) from good neighbours\n', nnz(drop));
end

fig = jneFigure(st.WidthSingle, st.WidthSingle*0.95, st);
ax = axes(fig, 'Position', [0.05 0.06 0.72 0.82]);
jneGridImage(ax, vals, R.channelMap, cmap, cl, st, ...
    'UpN',240, 'Mode',opt.Mode, 'NLevels',opt.NLevels, 'Interp',opt.Interp, ...
    'SmoothSigma',opt.SmoothSigma, 'ShowContacts',opt.ShowContacts);
title(ax, titleStr);
jneApplyAxes(ax, st);
% Vector colorbar (renders in Illustrator; MATLAB's colorbar exports as a raster pattern)
jneVectorColorbar(fig, [0.80 0.18 0.03 0.58], cmap, cl, opt.NLevels, cbLabel, st);

tag = char(opt.SaveTag); if isempty(tag), tag = defTag; end
if opt.Save, jneSave(fig, opt.SaveDir, tag, st); end
end


% =====================================================================
function [R, lab, goodMask] = resolveSession(D, dateStr)
% Use the cached representative R if this date is one of them; otherwise load
% the session's flashERP result on the fly. goodMask is the QC good-channel mask
% (all-true for cached reps, whose qc isn't stored in the slim cache).
if isfield(D,'flash') && ~isempty(D.flash)
    hit = find(string({D.flash.date}) == dateStr, 1);
    if ~isempty(hit)
        R = D.flash(hit).R; lab = D.flash(hit).label; goodMask = true(R.nChan,1); return;
    end
end
% release: the other Fig 7 sessions are in D.other (jneData), with their QC good-channel mask
hit = find(string({D.other.date}) == dateStr, 1);
if isempty(hit), error('fig7c_polarity_map:noSession','No Fig 7 session for %s', dateStr); end
R = D.other(hit).R; goodMask = R.goodMask;
lab = ageLabelFromDate(dateStr, D);
end

function lab = ageLabelFromDate(dateStr, D)
lab = char(dateStr);
if ~isfield(D,'implant') || isempty(D.implant), return; end
s = char(dateStr);
if numel(s)~=6 || any(~isstrprop(s,'digit')), return; end
dt = datetime(2000+str2double(s(1:2)), str2double(s(3:4)), str2double(s(5:6)));
ad = days(dt - D.implant);
if ~isfinite(ad), return; end
if ad < 56, lab = sprintf('%d wk', max(0,round(ad/7)));
else,        lab = sprintf('%d mo', round(ad/30.4)); end
end
