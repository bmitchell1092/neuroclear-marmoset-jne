function fig = fig7d_grid_waveforms(D, st, varargin)
% FIG 7d (JNE) — Trial-averaged evoked waveforms laid on the array grid at
% 2 wk / 2 mo / 4 mo on a SINGLE SHARED amplitude scale (honest amplitude
% comparison, consistent with the p2p heatmap). One representative contact is
% ringed in the grid; its mean +/- SEM waveform is produced SEPARATELY by
% jne_contactInset (one small figure per age) so it can be scaled/placed by hand.
%
% Pick the ringed / inset contact with 'InsetContact':
%   fig7d_grid_waveforms(D, st, 'InsetContact', 54)      % erp-row index (1..64)
%   fig7d_grid_waveforms(D, st, 'InsetContact', [8 4])   % [gridRow gridCol]
%   fig7d_grid_waveforms(D, st, 'InsetContact', [])      % [] = auto-pick strongest
if nargin < 1 || isempty(D),  D  = jneData();  end
if nargin < 2 || isempty(st), st = jneStyle(); end
p = inputParser;
p.addParameter('XLimMs', [-100 300], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('InsetContact', 8, @(x)isempty(x)||isscalar(x)||numel(x)==2); % [row col], erp-row, or []=auto
p.addParameter('SaveTag', 'fig7d_grid_waveforms', @(s)ischar(s)||isstring(s));
p.addParameter('SaveDir', jneOutDir(), @(s)ischar(s)||isstring(s));
p.addParameter('Save', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt = p.Results;

n = numel(D.flash);
chMap = D.flash(1).R.channelMap;

% ---- Shared vertical scale across panels ----
ampScale = 0;
for i = 1:n
    R = D.flash(i).R; sel = R.responsive;
    if ~any(sel), sel = true(R.nChan,1); end
    ampScale = max(ampScale, prctile(R.ampP2P(sel), 95));
end
if ~isfinite(ampScale) || ampScale<=0, ampScale = 100; end

% ---- Resolve ringed contact -> erp row + grid location ----
featCh = jneResolveContact(opt.InsetContact, D, chMap, n);
[fr, fc] = find(chMap == featCh-1, 1);
fprintf('[fig7d] ringed contact: erp-row %d  at grid (row %d, col %d)\n', featCh, fr, fc);

% ---- Grid panels (shared scale) ----
fig = jneFigure(st.WidthDouble, st.WidthDouble/n*1.28, st);
tl = tiledlayout(fig, 1, n, 'Padding','compact','TileSpacing','compact');
for i = 1:n
    ax = nexttile(tl); R = D.flash(i).R; [nr,~] = size(R.channelMap);
    jneWaveGrid(ax, R.erp, R.tMs, R.channelMap, ampScale, st.Black, opt.XLimMs, st);
    if ~isempty(fr)
        line(ax, fc, nr-fr+1, 'Marker','o', 'MarkerSize',9, 'LineWidth',1.2, ...
            'Color',st.Blue, 'MarkerFaceColor','none');
    end
    title(ax, D.flash(i).label, 'FontName',st.FontName, 'FontSize',st.FontSizeTitle, 'Interpreter','none');
end

if opt.Save, jneSave(fig, opt.SaveDir, char(opt.SaveTag), st); end
end
