function P = jneSpontEpochData(varargin)
% jneSpontEpochData  The Fig 8a spontaneous exemplar (and Supplementary Video 1).
%
% One eyes-open rest window (260710, window 148), re-cut with 500 ms padding either side so
% the band-pass edge transients fall outside the shown +/-250 ms. Per ms of the window:
%   Z      LFP 5-40 Hz (the detector's band), z-scored per contact
%   ph     generalized phase, 5-40 Hz (waveAnalyticFrames)
%   E      detected wave events (waveEvents, rho >= 0.37)
q = inputParser;
q.addParameter('RhoThr', 0.37);
q.parse(varargin{:}); o = q.Results;
dd = fullfile(jneDataDir(), 'fig8');
C = load(fullfile(dd, 'lightRest_exemplar.mat'));
gm = C.gridMask; tp = C.time_ms(:).'; half = 250;
core = tp >= -half & tp <= half; t = tp(core);
% the padded window's core is exactly that window of the session cube
W = load(fullfile(dd, sprintf('lightRest_%s.mat', extractBefore(C.session,7))), 'cube');
dv = max(abs(squeeze(C.cube(1,:,:,core)) - squeeze(W.cube(C.win,:,:,:))), [], 'all', 'omitnan');
assert(dv < 1e-6, 'jneSpontEpochData: exemplar does not match window %d of the session', C.win);
raw = squeeze(double(C.cube(1,:,:,:)));
for j = 1:size(raw,3)                                   % same hole fill as waveEpochMovie
    f = raw(:,:,j); if any(isnan(f(:))), raw(:,:,j) = fillDiscFrame(f); end
end
Z = zscore_independent(bandpass_filter(raw, 5, 40, 3, 1000)); Z = Z(:,:,core);
X = waveAnalyticFrames(squeeze(double(C.cube(1,:,:,:))), gm);
U = propagationField(X, gm, 1000);
[aD, aC] = waveMeanArrow(X, U, gm, find(core));
ph = angle(X(:,:,core));
Z(repmat(~gm,1,1,numel(t))) = NaN; ph(repmat(~gm,1,1,numel(t))) = NaN;
E = waveEvents(C.cube, gm, tp, 'RhoThr',o.RhoThr, 'detect_win_ms',[-half half], 'Verbose',false);
P = struct('session',C.session, 'win',C.win, 't',t, 'tElapsed',t - t(1), 'Z',single(Z), ...
    'ph',single(ph), 'arrowDir',aD, 'arrowCons',aC, 'E',E, 'gm',gm, 'pitchMm',C.pitch_mm, ...
    'cubePadded',single(C.cube), 'tPadded',tp);
end
