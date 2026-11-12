function P = jneFlashEpochData(varargin)
% jneFlashEpochData  The Fig 8a flash-evoked exemplar (and Supplementary Video 2).
%
% One trial (default 260724, trial 100), -300..400 ms around the flash. Same fields as
% jneSpontEpochData. Events are detected on the whole session cube, as in panel b, and
% this trial's are kept.
q = inputParser;
q.addParameter('Date', '260724');
q.addParameter('Epoch', 100);
q.addParameter('RhoThr', 0.37);
q.parse(varargin{:}); o = q.Results;
S = load(fullfile(jneDataDir(), 'fig8', sprintf('flash_%s.mat', o.Date)));
gm = S.gridMask; t = S.time_ms(:).';
raw = squeeze(double(S.cube(o.Epoch,:,:,:)));
for j = 1:size(raw,3)                                   % same hole fill as waveEpochMovie
    f = raw(:,:,j); if any(isnan(f(:))), raw(:,:,j) = fillDiscFrame(f); end
end
Z = zscore_independent(bandpass_filter(raw, 5, 40, 3, 1000));
X = waveAnalyticFrames(squeeze(double(S.cube(o.Epoch,:,:,:))), gm);
U = propagationField(X, gm, 1000);
[aD, aC] = waveMeanArrow(X, U, gm, 1:numel(t));
ph = angle(X);
Z(repmat(~gm,1,1,numel(t))) = NaN; ph(repmat(~gm,1,1,numel(t))) = NaN;
E = waveEvents(S.cube, gm, t, 'RhoThr',o.RhoThr, 'detect_win_ms',[-250 350], 'Verbose',false);
E = E(E.epoch == o.Epoch, :);
P = struct('session',S.session, 'epoch',o.Epoch, 't',t, 'Z',single(Z), 'ph',single(ph), ...
    'arrowDir',aD, 'arrowCons',aC, 'E',E, 'gm',gm, 'pitchMm',S.pitch_mm, ...
    'cube',single(S.cube(o.Epoch,:,:,:)));
end
