function fig = fig7e_p2p_progression(D, st, varargin)
% FIG 7e (JNE) — Peak-to-peak transient amplitude maps at 2 wk / 2 mo / 4 mo.
% Default 'Norm'='perpanel': each map is scaled to its own amplitude (peak µV
% labelled per panel) so the SHARED SPATIAL PATTERN is comparable across ages
% despite the amplitude gap (a single shared scale washes 2 wk out). 'Norm'='shared'
% uses one absolute colorbar ('CLim'). Contour mode (vector) is Illustrator-safe.
if nargin < 1 || isempty(D),  D  = jneData();  end
if nargin < 2 || isempty(st), st = jneStyle(); end
p = inputParser;
p.addParameter('TransientMs', [50 100], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Norm', 'perpanel', @(s)ischar(s)||isstring(s)); % 'perpanel' | 'shared'
p.addParameter('NormPct', 98, @isscalar);                      % per-panel scale percentile
p.addParameter('CLim', [0 65], @(x)isempty(x)||numel(x)==2);   % shared-mode colorbar (µV); []=auto
p.addParameter('Mode', 'contour', @(s)ischar(s)||isstring(s)); % 'contour' = vector (Illustrator-safe); 'smooth' = raster
p.addParameter('NLevels', 8, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ShowContacts', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('ShowPeakLabel', true, @(x)islogical(x)&&isscalar(x)); % per-panel peak µV (perpanel mode)
p.addParameter('InterpBadChans', false, @(x)islogical(x)&&isscalar(x)); % fill QC-bad contacts from good neighbours
p.addParameter('SaveTag', 'fig7e_p2p_progression', @(s)ischar(s)||isstring(s));
p.addParameter('SaveDir', jneOutDir(), @(s)ischar(s)||isstring(s));
p.addParameter('Save', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt = p.Results;

n = numel(D.flash);
P = cell(1,n); allv = [];
for i = 1:n
    R = D.flash(i).R;
    P{i} = jneWindowMetrics(R.erp, R.tMs, opt.TransientMs);
    if opt.InterpBadChans && isfield(R,'goodMask') && ~isempty(R.goodMask)
        drop = ~R.goodMask(:);
        if any(drop)
            P{i} = jneFillContacts(P{i}, R.channelMap, drop);
            fprintf('[fig7e] %s: filled %d contact(s)\n', char(D.flash(i).label), nnz(drop));
        end
    end
    allv = [allv; P{i}(R.responsive)]; %#ok<AGROW>
end

perpanel = strcmpi(opt.Norm, 'perpanel');
V = cell(1,n); scaleUV = nan(1,n);
if perpanel
    for i = 1:n
        sel = D.flash(i).R.responsive; if ~any(sel), sel = isfinite(P{i}); end
        s = prctile(P{i}(sel), opt.NormPct);
        if ~isfinite(s) || s <= 0, s = max(P{i}(isfinite(P{i}))); end
        scaleUV(i) = s;
        V{i} = P{i} / s;                         % normalised to its own scale
    end
    cl = [0 1]; cbLabel = 'Peak-to-peak (norm.)';
else
    V = P;
    if isempty(opt.CLim), cl = [0 prctile(allv, 99)]; else, cl = opt.CLim; end
    cbLabel = ['Peak-to-peak (' st.muV ')'];
end

fig = jneFigure(st.WidthDouble, st.WidthDouble/n*1.1, st);
tl = tiledlayout(fig, 1, n, 'Padding','compact','TileSpacing','compact');
try tl.OuterPosition = [0 0 0.90 1]; catch, end    % leave a right margin for the vector colorbar
for i = 1:n
    ax = nexttile(tl);
    jneGridImage(ax, V{i}, D.flash(i).R.channelMap, st.CMap_p2p, cl, st, ...
        'UpN',240, 'Mode',opt.Mode, 'NLevels',opt.NLevels, 'ShowContacts',opt.ShowContacts);
    title(ax, D.flash(i).label);
    jneApplyAxes(ax, st);
    if perpanel && opt.ShowPeakLabel
        text(ax, 0.03, 0.97, sprintf('%.0f %s', scaleUV(i), st.muV), 'Units','normalized', ...
            'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', ...
            'FontName',st.FontName, 'FontSize',st.FontSizeAnnot);
    end
end
% Vector colorbar (renders in Illustrator; MATLAB's colorbar exports as a raster pattern)
cbPos = [0.915 0.28 0.018 0.44];
if perpanel
    jneVectorColorbar(fig, cbPos, st.CMap_p2p, [0 1], opt.NLevels, cbLabel, st, 'Ticks',[0 0.5 1]);
else
    jneVectorColorbar(fig, cbPos, st.CMap_p2p, cl, opt.NLevels, cbLabel, st);
end

if opt.Save, jneSave(fig, opt.SaveDir, char(opt.SaveTag), st); end
end
