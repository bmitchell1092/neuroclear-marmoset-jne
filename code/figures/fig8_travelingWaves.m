function fig = fig8_travelingWaves(varargin)
% fig8_travelingWaves — Fig 8: traveling waves on the marmoset Neuro-Clear array.
%
% The NHP counterpart of the mouse Fig 6, built for two claims only: waves visibly cross the
% array, and they travel at physiological speeds.
%   a  one spontaneous wave (eyes-open rest) and one single-trial flash-evoked wave, each as
%      LFP (5-40 Hz, z-scored -- the SAME band the phase is computed in) over generalized
%      phase, per contact (no interpolation), six frames 10 ms apart. The arrow is the
%      DETECTED wave (waveEpochMovie 'phase+event'): its direction, length ~ rho, faded over
%      one cycle around the detection -- so an arrow appears only where the detector fired.
%      These are the frames of Supplementary Videos 1 and 2 (jne_suppVideos), same data.
%   b  detected wave events per trial around a full-field flash (bars, mean over sessions)
%      with the trial-averaged evoked LFP (red); bar = cluster-corrected rise above a
%      circular-shift null, session as the unit
%   c  per-event speed, rest vs flash-evoked (110-130 ms) waves (medians in the labels),
%      shaded = range reported in previous studies
%
% Data: jneSpontEpochData + jneFlashEpochData (a), jneFlashWaveData (b), jneSpeedData (c).
% Exemplars chosen for legibility from detected events -- the legend must say so.
%
% ⚠️ The in vivo band (0.1-0.8 m/s) must be checked against the cited papers.
%
% Name/value: 'Frames' -20:10:30 (spont, ms around the detection)
%             'FlashEvent' {'260724_Biru_fullFieldFlash_01', 100}  'FlashFrames' -5:5:20 (ms around
%             the flash detection: late-era evoked waves cross the array in ~8 ms, so 5 ms steps)
%             'EraFrom' '260622' (flash sessions used in b and c)  'ZLim' 3  'CycleMs' 70
%             'Band' [0.1 0.8]  'XLim' [-200 350]
%             'SaveDir' output/  'SaveTag' 'fig8_travelingWaves'  'Save' true
q = inputParser;
q.addParameter('Frames', -20:10:30);
% Fig 8 uses June 2026 onward (the early sessions are Fig 7's longevity span, with noisier
% wave-band phase). 260724 trial 100: best visible crest sweep of the late-era candidates
% (rho 0.75, crossing-time plane fit R2 0.80) -- Blake picked it over the day-10 trial.
q.addParameter('FlashEvent', {'260724_Biru_fullFieldFlash_01', 100});   % {session, trial} for jneFlashEpochData
q.addParameter('FlashFrames', -5:5:20);
q.addParameter('EraFrom', '260622');
q.addParameter('ZLim', 3);
q.addParameter('CycleMs', 70);
q.addParameter('Band', [0.1 0.8]);
q.addParameter('XLim', [-200 350]);
q.addParameter('SaveDir', jneOutDir());
q.addParameter('SaveTag', 'fig8_travelingWaves');
q.addParameter('Save', true);
q.parse(varargin{:}); o = q.Results;
st = jneStyle();
Ps = jneSpontEpochData(); F = jneFlashWaveData(); S = jneSpeedData();
Pf = jneFlashEpochData('Date',o.FlashEvent{1}(1:6), 'Epoch',o.FlashEvent{2});
% the exemplar event in each: the strongest detection (spont: +13 ms, rho 0.84 -- the
% demo window's other detections are weaker; flash: +123 ms, rho 0.75, at the event-rate peak)
[~, i1] = max(Ps.E.rho); [~, i2] = max(Pf.E.rho);
evS = Ps.E(i1,:); evF = Pf.E(i2,:);

W = st.WidthDouble; H = 3.7; Hd = 0;                    % (wave-DAQ's version has an optional longevity row)
fig = jneFigure(W, H + Hd, st);
hsvMap = hsv(256); ARROW_MAX = 5.0; red = [0.80 0.10 0.10];
[nr, nc] = size(Ps.gm); nf = numel(o.Frames); assert(numel(o.FlashFrames) == nf, 'fig8: both rows need the same frame count');
fw = 0.088; fg = 0.006; x0 = 0.075; fh = fw*W/H * nr/nc;     % frame size keeps the array's aspect
rowY = [0.745 0.565 0.27 0.09];                                % spont LFP, spont phase, flash LFP, flash phase

% ---- a  LFP over phase, per contact ----
grp = {Ps, evS, 'Spontaneous', 'time from wave detection'; ...
       Pf, evF, 'Flash-evoked', 'time from flash'};
for g = 1:2
    P = grp{g,1}; ev = grp{g,2}; fr = o.Frames; if g == 2, fr = o.FlashFrames; end
    for f = 1:nf
        tk = ev.time_ms + fr(f); k = find(P.t >= tk, 1);
        x = x0 + (f-1)*(fw+fg);
        % LFP row
        ax = axes(fig, 'Position', [x, rowY(2*g-1), fw, fh]);
        cellPatches(ax, double(P.Z(:,:,k)), P.gm, st.CMap_div, [-o.ZLim o.ZLim]);
        if g == 1, lbl = sprintf('%+d ms', fr(f)); else, lbl = sprintf('%d ms', round(tk)); end
        title(ax, lbl, 'FontSize',st.FontSizeAnnot, 'FontWeight','normal', 'Interpreter','none');
        if f == 1, rowLabel(ax, 'LFP', nr, st); end
        % phase row + the detected-wave arrow
        ax = axes(fig, 'Position', [x, rowY(2*g), fw, fh]);
        cellPatches(ax, double(P.ph(:,:,k)), P.gm, hsvMap, [-pi pi]);
        [u, v, w] = eventArrow(P.E, P.t(k), o.CycleMs, ARROW_MAX);
        if w > 0.02
            quiver(ax, (nc+1)/2 - u/2, (nr+1)/2 - v/2, u, v, 0, 'Color',[0 0 0 max(w,0.05)], ...
                'LineWidth',0.8 + 1.4*w, 'MaxHeadSize',1.2);
        end
        if f == 1, rowLabel(ax, 'Phase', nr, st); end
        if g == 2 && f == nf                            % 1 mm scale bar
            plot(ax, [nc+0.5-1/P.pitchMm, nc+0.5], [nr+0.9 nr+0.9], 'k-', 'LineWidth',1.5, 'Clipping','off');
            text(ax, nc+0.5-0.5/P.pitchMm, nr+1.1, '1 mm', 'HorizontalAlignment','center', ...
                'VerticalAlignment','top', 'FontSize',st.FontSizeAnnot-1, 'Clipping','off');
        end
    end
    wa = nf*(fw+fg) - fg;                                % group header over its two rows
    annotation(fig, 'textbox', [x0 rowY(2*g-1)+fh+0.045 wa 0.05], 'String', ...
        sprintf('%s  (%s)', grp{g,3}, grp{g,4}), 'EdgeColor','none', 'FontSize',st.FontSizeAnnot, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end
xr = x0 + nf*(fw+fg);
jneVectorColorbar(fig, [xr, rowY(1)+0.1*fh, 0.008, 0.8*fh], st.CMap_div, [-o.ZLim o.ZLim], 12, '', st, ...
    'Ticks',[-o.ZLim 0 o.ZLim], 'TickLabels',{sprintf('%g',-o.ZLim), '0', sprintf('%g',o.ZLim)});
jneVectorColorbar(fig, [xr, rowY(2)+0.1*fh, 0.008, 0.8*fh], hsvMap, [-pi pi], 24, '', st, ...
    'Ticks',[-pi 0 pi], 'TickLabels',{['-' char(960)], '0', char(960)});

% ---- b  wave events per trial around the flash ----
R = F.rate; xb = 0.755; wb = 0.185;
ax = axes(fig, 'Position', [xb 0.62 wb 0.29]); hold(ax,'on');
yyaxis(ax,'left');
bar(ax, R.t, R.mean, 1, 'FaceColor',[0.72 0.72 0.72], 'EdgeColor','none');
yl = [0 max(R.mean)*1.25]; ylim(ax, yl);
Cl = R.clusters;
for i = 1:height(Cl)
    if Cl.p(i) < 0.05 && Cl.sign(i) > 0
        plot(ax, [Cl.t1(i) Cl.t2(i)], yl(2)*0.95*[1 1], '-', 'Color','k', 'LineWidth',2);
    end
end
ylabel(ax, 'events / trial'); ax.YColor = [0.25 0.25 0.25];
yyaxis(ax,'right');
erp = squeeze(mean(F.erp.perSession, [1 2 4], 'omitnan'));
plot(ax, F.erp.t, erp, '-', 'Color',red, 'LineWidth',st.LineWidthData);
ylabel(ax, ['LFP (' st.muV ')']); ax.YColor = red;
xlim(ax, o.XLim); xline(ax, 0, ':', 'Color',st.Black, 'LineWidth',st.LineWidthThin);
xlabel(ax, 'time from flash (ms)');
title(ax, sprintf('%d sessions, %d trials', R.nSess, R.nTrials), 'FontSize',st.FontSizeAnnot, ...
    'FontWeight','normal', 'Interpreter','none');
jneApplyAxes(ax, st);

% ---- c  speed distributions ----
% flash-evoked speeds from the same era as b (all lightRest days in the LUT are Jul-Aug 2026)
E = S.E; xl = [0.02 3];
E = E(E.surface ~= "flash" | double(extractBefore(E.session,7)) >= str2double(o.EraFrom), :);
ax = axes(fig, 'Position', [xb 0.11 wb 0.29]); hold(ax,'on');
xi = linspace(log10(xl(1)), log10(xl(2)), 300);
sg = {"spont","spont",[0.15 0.15 0.15],'Spontaneous'; "flash","evoked",[0.85 0.33 0.10],'Flash-evoked'};
dens = cell(1,2); pk = 0;
for g = 1:2
    v = E.speed_m_s(E.surface == sg{g,1} & E.window == sg{g,2});
    dens{g} = ksdensity(log10(v), xi); pk = max(pk, max(dens{g}));
end
patch(ax, o.Band([1 2 2 1]), [0 0 1.6 1.6]*pk, [0.88 0.88 0.88], 'EdgeColor','none');
for g = 1:2
    v = E.speed_m_s(E.surface == sg{g,1} & E.window == sg{g,2});
    plot(ax, 10.^xi, dens{g}, '-', 'Color',sg{g,3}, 'LineWidth',st.LineWidthData);
    text(ax, 0.022, (1.68 - 0.15*g)*pk, sprintf('%s  %.2f m/s', sg{g,4}, median(v)), ...     % empty slow tail
        'Color',sg{g,3}, 'FontSize',st.FontSizeAnnot-1, 'VerticalAlignment','top', 'Interpreter','none');
end
text(ax, sqrt(prod(o.Band)), 0.03*pk, 'reported', 'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', 'FontSize',st.FontSizeAnnot-1, 'Color',[0.35 0.35 0.35]);
set(ax, 'XScale','log', 'XLim',xl, 'YLim',[0 1.6*pk], 'YTick',[], 'XTick',[0.03 0.1 0.3 1 3], 'XTickLabelRotation',0);
xlabel(ax, 'wave speed (m/s)'); ylabel(ax, 'density');
jneApplyAxes(ax, st);

% ---- panel letters (a-c) ----
lt = {'a', [0.005 0.93]; 'b', [xb-0.075 0.93]; 'c', [xb-0.075 0.42]};
for i = 1:3
    annotation(fig, 'textbox', [lt{i,2} 0.04 0.06], 'String',lt{i,1}, 'EdgeColor','none', ...
        'FontSize',st.FontSizeTitle+1, 'FontWeight','bold', 'VerticalAlignment','bottom');
end

if o.Save, jneSave(fig, o.SaveDir, o.SaveTag, st); end
end

% =====================================================================
function cellPatches(ax, M, gm, cmap, cl)
% one flat square per real contact (vector; no interpolation) -- as the video draws them
[nr,nc] = size(gm); V = []; Fc = []; Cd = [];
for rr = 1:nr
    for cc = 1:nc
        if ~gm(rr,cc) || ~isfinite(M(rr,cc)), continue; end
        n0 = size(V,1);
        V = [V; cc-.5 rr-.5; cc+.5 rr-.5; cc+.5 rr+.5; cc-.5 rr+.5]; %#ok<AGROW>
        Fc = [Fc; n0+(1:4)]; %#ok<AGROW>
        x = min(max((M(rr,cc)-cl(1))/diff(cl), 0), 1);
        Cd = [Cd; cmap(1 + round(x*(size(cmap,1)-1)), :)]; %#ok<AGROW>
    end
end
patch(ax, 'Faces',Fc, 'Vertices',V, 'FaceVertexCData',Cd, 'FaceColor','flat', 'EdgeColor','none');
hold(ax,'on'); axis(ax,'image'); set(ax, 'YDir','reverse', 'XTick',[], 'YTick',[], ...
    'XLim',[0.5 nc+0.5], 'YLim',[0.5 nr+0.5], 'XColor','none', 'YColor','none');
end

function [u, v, w] = eventArrow(Ev, tk, cycleMs, arrowMax)
% waveEpochMovie 'phase+event', verbatim: the detected wave's direction (detector frame,
% top-CW -> screen), length ~ rho, raised-cosine weight over one cycle around the detection
u = 0; v = 0; w = 0;
for e = 1:height(Ev)
    dt = tk - Ev.time_ms(e);
    if abs(dt) > cycleMs/2, continue; end
    we = 0.5*(1 + cos(pi*dt/(cycleMs/2)));
    if we <= w, continue; end
    g = exp(1i*deg2rad(90 - Ev.dirDeg(e))); d = [real(g), -imag(g)]; d = d/norm(d);
    L = Ev.rho(e) * arrowMax; u = d(1)*L; v = d(2)*L; w = we;
end
end

function rowLabel(ax, s, nr, st)
text(ax, -0.2, (nr+1)/2, s, 'Rotation',90, 'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', 'FontSize',st.FontSizeAnnot, 'Interpreter','none');
end
