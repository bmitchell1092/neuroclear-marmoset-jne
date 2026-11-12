function figs = jne_contactInset(D, st, varargin)
% jne_contactInset — standalone SEM-inset figures (one per representative session)
% for a SINGLE contact: trial-averaged flash VEP, mean +/- SEM, peak-to-peak
% normalised within the box. Companion to fig7d_grid_waveforms (panel d):
% fig7d rings the contact on the array; these small figures are exported
% separately so you can scale and place them by hand (e.g. in Illustrator).
% Aspect ~ the fig7d inset (slightly wider than tall).
%
% Pick the contact with 'InsetContact' (same forms as fig7d):
%   jne_contactInset(D, st, 'InsetContact', [8 4])   % [gridRow gridCol]
%   jne_contactInset(D, st, 'InsetContact', 54)      % erp-row index
%   jne_contactInset(D, st, 'InsetContact', [])      % auto-pick strongest
%
% One PNG+SVG per age is written as <SaveTag>_<age>_r<row>c<col>.* (SaveTag
% default 'fig7d_inset'). 'Width'/'Height' set the figure size in inches.
if nargin < 1 || isempty(D),  D  = jneData();  end
if nargin < 2 || isempty(st), st = jneStyle(); end
p = inputParser;
p.addParameter('InsetContact', 8, @(x)isempty(x)||isscalar(x)||numel(x)==2);
p.addParameter('InsetXLimMs', [-50 250], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Width', 1.8, @isscalar);          % inches (aspect ~1.2, like the old inset)
p.addParameter('Height', 1.5, @isscalar);
p.addParameter('ShowAge', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('SaveTag', 'fig7d_inset', @(s)ischar(s)||isstring(s));
p.addParameter('SaveDir', jneOutDir(), @(s)ischar(s)||isstring(s));
p.addParameter('Save', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt = p.Results;

n = numel(D.flash);
chMap = D.flash(1).R.channelMap;
featCh = jneResolveContact(opt.InsetContact, D, chMap, n);
[fr, fc] = find(chMap == featCh-1, 1);
fprintf('[jne_contactInset] contact erp-row %d at grid (row %d, col %d)\n', featCh, fr, fc);

figs = gobjects(1,n);
for i = 1:n
    fig = jneFigure(opt.Width, opt.Height, st);
    ax  = axes('Parent',fig, 'Position',[0.05 0.05 0.90 0.90]);
    jneContactWave(ax, D.flash(i).R, featCh, opt.InsetXLimMs, st);
    if opt.ShowAge
        text(ax, 0.96, 0.96, char(D.flash(i).label), 'Units','normalized', ...
            'HorizontalAlignment','right', 'VerticalAlignment','top', 'Interpreter','none', ...
            'FontName',st.FontName, 'FontSize',st.FontSizeAnnot-1);
    end
    figs(i) = fig;
    if opt.Save
        age = regexprep(char(D.flash(i).label), '\s', '');
        tag = sprintf('%s_%s_r%dc%d', char(opt.SaveTag), age, fr, fc);
        jneSave(fig, opt.SaveDir, tag, st);
    end
end
end
