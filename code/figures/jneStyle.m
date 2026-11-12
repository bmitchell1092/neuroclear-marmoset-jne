function st = jneStyle()
% jneStyle
%
% THE uniform style for the JNE paper figure set (Dr. Kazum manuscript). Single
% source of truth — edit values here and every figure that uses jneFigure /
% jneApplyAxes / the returned struct updates. Consumed across analysis/JNE_paper/.
%
% Returns a struct of fonts, line widths, colors, colormaps, panel sizes, and
% export settings. Keep this the ONLY place these are defined.

st = struct();

% ---- Fonts (pt) ----
st.FontName      = 'Arial';
st.FontSizeTick  = 8;      % axis tick labels
st.FontSizeLabel = 9;      % x/y labels
st.FontSizeTitle = 10;     % panel titles
st.FontSizeAnnot = 8;      % in-axes text / legends

% ---- Lines / axes ----
st.LineWidthData = 1.2;    % primary data traces
st.LineWidthThin = 0.5;    % secondary / many-channel traces
st.LineWidthAxis = 0.75;   % axis rule + tick width
st.Box           = 'on';     % full frame around axes (collaborator preference)
st.TickDir       = 'in';     % ticks inward so they don't collide with the box
st.TickLength    = [0.02 0.02];
st.GridAlpha     = 0.12;

% ---- Colors ----
st.Blue = [0.10 0.45 0.80];
st.Gray = [0.70 0.70 0.70];
st.Black = [0 0 0];
% Age palette (2 wk / 2 mo / 4 mo) — pastel blue / green / orange, even brightness
st.AgeColors = [0.45 0.68 0.86;
                0.52 0.78 0.55;
                0.97 0.74 0.42];
st.AgeLabels = ["2 wk" "2 mo" "4 mo"];

% ---- Unicode glyphs (used with Interpreter 'none' so SVG text stays editable) ----
% TeX markup (\muV, \pm, ^2) makes MATLAB split labels into per-glyph runs in SVG,
% which Illustrator can't edit as whole lines. Use these literals + 'none' instead.
st.muV       = [char(181) 'V'];                          % µV
st.pm        = char(177);                                % ±
st.psdYLabel = ['PSD (dB ' char(181) 'V' char(178) '/Hz)']; % PSD (dB µV²/Hz)

% ---- Colormaps ----
st.CMap_p2p     = turbo(256);    % peak-to-peak magnitude maps (figs 3,4)
st.CMap_latency = parula(256);   % peak-latency maps (fig 5)
tdiv = linspace(0,1,256).';      % diverging blue-white-red for signed voltage
st.CMap_div     = interp1([0 .5 1], [0.13 0.40 0.67; 1 1 1; 0.70 0.09 0.17], tdiv);
st.NaNColor     = [0.9 0.9 0.9]; % absent/masked grid cells

% ---- Panel geometry (inches) ----
st.WidthSingle = 3.5;   % single-column
st.WidthDouble = 7.2;   % double-column
st.RowHeight   = 2.4;   % default panel-row height

% ---- Export ----
st.ExportDPI = 300;
end
