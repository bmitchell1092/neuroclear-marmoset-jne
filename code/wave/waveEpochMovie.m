function waveEpochMovie(cube, gridMask, time_ms, sIdx, outFile, varargin)
% waveEpochMovie   Render ONE epoch's traveling wave as an mp4 -- any trigger, or none.
%
% Two panels over time: z-scored broadband LFP (left) and generalized phase (right).
% For a non-rectangular array the holes are filled for the GP/gradient but masked
% (transparent) in display, so you only ever see real electrodes.
%
% ⚠ DEFAULT STYLE CHANGED 2026-09-08 to 'phase+mean' -- GP phase plus ONE centred arrow.
% The previous default drew one arrow per electrode ('phase+field'); readers consistently
% found 64 arrows too busy to follow, and most of them carry no detection anyway. Pass
% 'Style','phase+field' for the old look. Movies rendered before this date show the old
% style; nothing about the DATA changed, only the overlay.
%
% INPUTS
% ------
% cube      nEpoch x nRow x nCol x T (NaN at non-electrode cells)
% gridMask  nRow x nCol logical (true = real electrode)
% time_ms   1 x T (ms, 0 = trigger)
% sIdx      epoch index to render
% outFile   output .mp4 path
% Name/value: 'Fs'1000 'gp_lp'0 'gp_band'[5 40] 'z_band'[0.4 40] 'FPS'12
%             'WinMs'[-100 250] 'ZLim'3 'Quiver'true
%             'Style' 'phase+mean' | 'phase' | 'phase+field' | 'phase+event' | 'residual'
%             'Events' [] (required by 'phase+event')  'CycleMs' 70  'ResidLim' 0.8
%   gp_band/gp_lp MUST match waveSourcePoints (5-40 Hz, lp 0) or the movie shows a
%   different phase field from the one the detections came from -- see the note below.
%             'SidePanel'[] 'SidePanelTitle'''
%
% Display is in natural channelMap orientation (row 1 = top, col 1 = left; NO
% transpose) -- confirm this matches the array's physical/anatomical orientation
% before reading propagation directions off the movie.
%
% SIDE PANEL (companion time series synchronised to the wave)
% ----------------------------------------------------------
% 'SidePanel' takes a function handle and adds a third tile that is redrawn every
% frame, so any co-recorded signal can play out beside the wave (screen-space
% gaze, pupil, a stimulus trace, a single-channel LFP). The engine stays agnostic
% about WHAT is drawn -- it supplies only the axis and the current frame time; the
% caller closes over its own data. Signature:
%
%   state = fn(ax, t_now_ms, state)      % state is [] on the FIRST call
%
% Draw static content (paths, targets, limits) on the first call and return the
% handles you need; on later calls just move the time-varying bits (cursor, dot).
% Returning [] every call is legal but redraws from scratch each frame.
% Example consumer: analysis/deliverables/saccadeTalk/slide5_waveGaze.m
%
% See also: waveGP, waveSourcePoints, accumulateSaccadeCube

p = struct('Fs',1000,'gp_lp',0,'gp_band',[5 40],'z_band',[0.4 40],'FPS',12,'WinMs',[-100 250],'ZLim',3,'Quiver',true, ...
           'Style','phase+mean', 'Events',[], 'CycleMs',70, 'ResidLim',0.8, ...  % right-panel rendering; see STYLE note
           'Trigger','', 'EpochLabel','epoch', 'TNote','', ...  % Trigger: '' = untriggered (time counts up); otherwise the event name shown at t=0
           'SidePanel',[],'SidePanelTitle','');
given = varargin(1:2:end);
for i=1:2:numel(varargin), p.(varargin{i})=varargin{i+1}; end

% BACKWARD COMPATIBILITY FOR 'Quiver'. The default Style changed to 'phase+mean' on
% 2026-09-08 (the 64-arrow field was judged unreadable). A caller that explicitly asked
% for 'Quiver',true wanted the per-electrode field, so honour that rather than silently
% handing it a different picture -- slide5_waveGaze does exactly this. An explicit
% 'Style' always wins; 'Quiver' only steers when Style was not given.
if ~any(strcmpi(given,'Style')) && any(strcmpi(given,'Quiver'))
    if p.Quiver, p.Style = 'phase+field'; else, p.Style = 'phase'; end
end
gm = logical(gridMask); [nRow,nCol] = size(gm);

raw = squeeze(cube(sIdx,:,:,:));                 % nRow x nCol x T (NaN corners)
T = size(raw,3);
for k=1:T, fr=raw(:,:,k); if any(isnan(fr(:))), raw(:,:,k)=fillDiscFrame(fr); end, end

% BANDPASS BEFORE GENERALIZED PHASE -- this must match what the DETECTOR does, and
% until 2026-08-31 it did not. The movie called generalized_phase on the raw cube with
% gp_lp = 0.5 and no bandpass, while waveSourcePoints bandpasses 5-40 Hz with
% gp_lp = 0. Measured on flash epoch 28: the movie's phase advanced at an implied
% **90.9 Hz** against the detector's **13.9 Hz** -- and since LFP_filt is already
% limited to 0.4-40 Hz, 91 Hz cannot be an oscillation at all; it was phase slipping
% erratically. The phase panel raced while the LFP panel (dominant 1.4 Hz) crawled,
% which is exactly what it looked like (John, 2026-08-31). Every movie rendered before
% this fix shows a phase field that is NOT the one the statistics were computed on.
xf = bandpass_filter(double(raw), p.gp_band(1), p.gp_band(2), 4, p.Fs);
X  = generalized_phase(xf, p.Fs, p.gp_lp);  Ph = angle(X);
% The LFP panel stays BROADBAND on purpose: it is there to show the evoked deflection,
% which lives below the wave band. The two panels therefore move at different rates by
% design -- both bands are printed in the panel titles so that is legible rather than
% confusing.
bp = bandpass_filter(double(raw), p.z_band(1), p.z_band(2), 3, p.Fs);
Z  = zscore_independent(bp);

% mask corners for display
alpha = double(gm);
frames = find(time_ms>=p.WinMs(1) & time_ms<=p.WinMs(2));
[Xg,Yg] = meshgrid(1:nCol,1:nRow);

% precompute the sign-corrected PROPAGATION field (not raw gradient) + ONE global
% arrow scale so lengths are comparable across frames and never rescale.
% Drawn as (real(U), -imag(U)) for the YDir='reverse' axis => arrows point where the
% wave travels (matches waveGP meanDirDeg / waveFieldSnapshot).
Uprop = propagationField(X, gm, p.Fs);
GX = nan(nRow,nCol,numel(frames)); GY = GX;
for i = 1:numel(frames)
    u = real(Uprop(:,:,frames(i))); v = -imag(Uprop(:,:,frames(i)));
    u(~gm)=NaN; v(~gm)=NaN; GX(:,:,i)=u; GY(:,:,i)=v;
end
mag = hypot(GX,GY); q95 = prctile(mag(~isnan(mag)),95);
qscale = 0.7 / max(q95, eps);            % typical arrow ~0.7 cells

% ---- RIGHT-PANEL STYLE ----------------------------------------------------------
% 'phase+field'  GP phase + one arrow per electrode  (the pre-2026-09-08 default)
% 'phase'        GP phase alone
% 'phase+mean'   GP phase + ONE centred arrow, every frame; length = directional
%                CONSISTENCY, so it shrinks when the array disagrees   <-- DEFAULT
% 'phase+event'  GP phase + ONE centred arrow drawn ONLY while a wave is DETECTED;
%                length proportional to that event's rho, raised-cosine fade over one
%                cycle. Needs 'Events' (the rows of the waveEvents table for THIS
%                epoch). This is the only style that distinguishes detected moments
%                from undetected ones, which is the distinction the analysis rests on.
% 'residual'     phase MINUS the array's own mean phase each frame. With only ~0.12 of
%                a wavelength across this array (median 49 mm vs ~9 mm of array), the
%                common-mode oscillation dominates and hides the spatial ramp; removing
%                it exposes the lead/lag structure. ⚠ It is a RELATIVE measure -- some
%                sites are red and some blue even in noise -- so it illustrates a wave,
%                it never demonstrates one. Keep its colourbar visible.
STYLE = lower(string(p.Style));
assert(ismember(STYLE, ["phase+field","phase","phase+mean","phase+event","residual"]), ...
    'waveEpochMovie: unknown Style ''%s''', p.Style);
isResid = STYLE == "residual";

% Display cube: phase, or phase with the per-frame spatial mean removed.
D = Ph;
if isResid
    for i = 1:numel(frames)
        k = frames(i); v = X(:,:,k); v = v(gm & isfinite(v));
        if isempty(v), continue; end
        f = angle(exp(1i*(Ph(:,:,k) - angle(mean(v))))); f(~gm) = NaN;
        D(:,:,k) = f;
    end
end

% ONE-ARROW styles, precomputed per frame: direction (unit), length, and fade weight.
cxA = (1+nCol)/2; cyA = (1+nRow)/2; ARROW_MAX = 5.0;    % cells at full length
[AU, AV, AW] = deal(zeros(1,numel(frames)));
if STYLE == "phase+mean"
    Amp = abs(X);
    for i = 1:numel(frames)
        k = frames(i); u = Uprop(:,:,k); a = Amp(:,:,k); u(~gm) = NaN; a(~gm) = NaN;
        gsum = sum(u(:).*a(:),'omitnan') / (sum(a(:),'omitnan') + eps);
        un = u(isfinite(u) & abs(u) > 0);
        if isempty(un) || ~isfinite(gsum) || abs(gsum) == 0, continue; end
        cons = abs(mean(un./abs(un)));                  % directional consistency 0-1
        d = [real(gsum), -imag(gsum)] / abs(gsum);      % -imag: YDir is reverse
        AU(i) = d(1)*cons*ARROW_MAX; AV(i) = d(2)*cons*ARROW_MAX; AW(i) = cons;
    end
elseif STYLE == "phase+event"
    Ev = p.Events;
    assert(~isempty(Ev) && istable(Ev) && all(ismember({'time_ms','rho','dirDeg'}, ...
        Ev.Properties.VariableNames)), ...
        ['waveEpochMovie: Style ''phase+event'' needs ''Events'' -- the waveEvents rows ' ...
         'for this epoch, with time_ms, rho and dirDeg.']);
    for i = 1:numel(frames)
        tk = time_ms(frames(i)); best = 0;
        for e = 1:height(Ev)
            dtE = tk - Ev.time_ms(e);
            if abs(dtE) > p.CycleMs/2, continue; end
            w = 0.5*(1 + cos(pi*dtE/(p.CycleMs/2)));    % raised cosine over one cycle
            if w <= best, continue; end
            best = w;
            g = exp(1i*deg2rad(90 - Ev.dirDeg(e)));      % top-CW -> math, as waveGP
            d = [real(g), -imag(g)]; d = d/norm(d);
            L = Ev.rho(e) * ARROW_MAX;
            AU(i) = d(1)*L; AV(i) = d(2)*L; AW(i) = w;
        end
    end
end

hasPanel = ~isempty(p.SidePanel);
figW = 1150 + 430*hasPanel;               % third tile widens the frame, same height
fig = figure('Color','w','Position',[60 60 figW 470],'Visible','off');
tl  = tiledlayout(fig,1,2+hasPanel,'Padding','compact','TileSpacing','compact');
axZ = nexttile(tl); hZ = imagesc(axZ, Z(:,:,frames(1)),'AlphaData',alpha);
colormap(axZ, redblue_local(256)); clim(axZ,[-p.ZLim p.ZLim]); colorbar(axZ,'southoutside');
axis(axZ,'image');
set(axZ,'Color',[.8 .8 .8],'XTick',[],'YTick',[], ...
        'XLim',[0.5 nCol+0.5],'YLim',[0.5 nRow+0.5],'XLimMode','manual','YLimMode','manual');
title(axZ, sprintf('z-LFP  (%.1f-%g Hz)', p.z_band(1), p.z_band(2)));
axP = nexttile(tl); hold(axP,'on');
hP = imagesc(axP, D(:,:,frames(1)),'AlphaData',alpha);
if isResid
    colormap(axP, redblue_local(256)); clim(axP,[-p.ResidLim p.ResidLim]);
    colorbar(axP,'southoutside');
else
    colormap(axP, hsv); clim(axP,[-pi pi]);
    colorbar(axP,'southoutside','Ticks',[-pi 0 pi],'TickLabels',{'-\pi','0','\pi'});
end
axis(axP,'image');
set(axP,'Color',[.8 .8 .8],'XTick',[],'YTick',[],'YDir','reverse', ...
        'XLim',[0.5 nCol+0.5],'YLim',[0.5 nRow+0.5],'XLimMode','manual','YLimMode','manual');
titles = struct('phase_field','GP phase + propagation  (%g-%g Hz, as detected)', ...
                'phase','GP phase  (%g-%g Hz, as detected)', ...
                'phase_mean','GP phase + mean propagation  (%g-%g Hz)   arrow length = consistency', ...
                'phase_event','GP phase + DETECTED wave  (%g-%g Hz)   arrow length ~ rho', ...
                'residual','phase MINUS array mean  (%g-%g Hz)   red = ahead, blue = behind');
title(axP, sprintf(titles.(strrep(char(STYLE),'+','_')), p.gp_band(1), p.gp_band(2)));

% Per-electrode field: kept for 'phase+field' only. 'Quiver',false still suppresses it,
% so existing callers behave exactly as before.
hQ = quiver(axP, Xg, Yg, qscale*GX(:,:,1), qscale*GY(:,:,1), 'k', 'AutoScale','off');
if ~p.Quiver || STYLE ~= "phase+field", set(hQ,'Visible','off'); end

% The single centred arrow, created once and updated per frame.
hA = quiver(axP, cxA, cyA, 0, 0, 0, 'Color',[0 0 0], 'LineWidth',3, 'MaxHeadSize',1.4, ...
            'AutoScale','off', 'Visible','off');

% Optional third tile: caller-supplied companion signal, redrawn per frame.
axS = []; panelState = [];
if hasPanel
    axS = nexttile(tl);
    if ~isempty(p.SidePanelTitle), title(axS, p.SidePanelTitle); end
end
st = sgtitle(fig,'');

v = VideoWriter(outFile,'MPEG-4'); v.FrameRate=p.FPS; v.Quality=100; open(v);
for i = 1:numel(frames)
    t = frames(i);
    hZ.CData = Z(:,:,t);
    hP.CData = D(:,:,t);
    if p.Quiver && STYLE == "phase+field"
        set(hQ,'UData',qscale*GX(:,:,i),'VData',qscale*GY(:,:,i));
    end
    if AW(i) > 0.02
        % Fade by alpha AND by line width: alpha alone is easy to miss on a projector,
        % and width alone cannot go to nothing. Arrow is centred, so it reads as
        % "the array is doing this" rather than claiming an origin -- the source point
        % is rim-biased on an array this much smaller than a wavelength.
        set(hA, 'UData',AU(i), 'VData',AV(i), 'Visible','on', ...
                'Color',[0 0 0 max(AW(i),0.05)], 'LineWidth',1.5 + 3.5*AW(i));
    else
        set(hA, 'Visible','off');
    end
    if hasPanel, panelState = p.SidePanel(axS, time_ms(t), panelState); end
    st.String = titleFor(p, sIdx, time_ms(t), time_ms(frames(1)));
    drawnow; writeVideo(v, getframe(fig));
end
close(v); close(fig);
fprintf('waveEpochMovie: wrote %s (%d frames)\n', outFile, numel(frames));
end

% ----------------------------------------------------------------------
function c = redblue_local(n)
if exist('customRedBlue','file')==2, c = customRedBlue(n); return; end
h = linspace(0,1,n)'; mid = round(n/2);
c = ones(n,3);
c(1:mid,1) = h(1:mid)/h(mid); c(1:mid,2) = h(1:mid)/h(mid);       % blue->white
c(mid+1:end,2) = flipud(h(1:n-mid))/h(mid); c(mid+1:end,3) = flipud(h(1:n-mid))/h(mid); % white->red
c = min(max(c,0),1);
end

% ----------------------------------------------------------------------
function s = titleFor(p, idx, tNow, tFirst)
% TRIGGERED vs UNTRIGGERED titling. With a trigger, t is signed and relative to the
% event ("+123 ms" after flash onset). Without one there is no zero to be relative
% TO, so showing "+/-" ms would invent an event that does not exist -- spontaneous
% epochs are windows, not responses. Time counts UP from the start of the window
% instead, which is what it actually is.
if isempty(p.Trigger)
    s = sprintf('%s %d   t = %.0f ms (elapsed; no trigger)', p.EpochLabel, idx, tNow - tFirst);
else
    s = sprintf('%s %d   t = %+d ms from %s onset', p.EpochLabel, idx, round(tNow), p.Trigger);
end
if ~isempty(p.TNote), s = sprintf('%s   %s', s, p.TNote); end
end
