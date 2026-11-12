function files = jne_suppVideos(varargin)
% jne_suppVideos — Supplementary Videos 1 and 2 for Fig 8: the full spontaneous fixation
% window and the full flash trial behind the Fig 8a frames, rendered by waveEpochMovie from
% the SAME data the frames use (jneSpontEpochData / jneFlashEpochData).
%
% Left: LFP 5-40 Hz, z-scored per contact (the detector's band, so it matches the phase).
% Right: generalized phase, with ONE arrow drawn only while a wave is DETECTED ('phase+event':
% direction of the detected wave, length proportional to rho, faded over one cycle).
% One video frame per ms at 12 fps, i.e. ~83x slower than real time -- the legend must say so.
%
% Name/value: 'SaveDir' output/  'FPS' 12
% See also: fig8_travelingWaves, waveEpochMovie
q = inputParser;
q.addParameter('SaveDir', jneOutDir());
q.addParameter('FPS', 12);
q.parse(varargin{:}); o = q.Results;
S = jneFlipEpoch(jneSpontEpochData()); F = jneFlipEpoch(jneFlashEpochData());   % figure orientation
files = {fullfile(o.SaveDir, 'SupplementaryVideo1_spontaneousWave.mp4'), ...
         fullfile(o.SaveDir, 'SupplementaryVideo2_flashEvokedWave.mp4')};
% waveEpochMovie titles each frame '<EpochLabel> <epoch index>'; the cube is passed as epoch
% 1 (video 1) or duplicated into epoch 2 (video 2) so the title reads "Supplementary Video N"
waveEpochMovie(double(S.cubePadded), S.gm, S.tPadded, 1, files{1}, 'z_band',[5 40], ...
    'Style','phase+event', 'Events',S.E, 'WinMs',[min(S.t) max(S.t)], 'FPS',o.FPS, ...
    'EpochLabel','Supplementary Video', 'TNote','spontaneous, eyes-open rest');
c2 = double(F.cube); c2 = cat(1, c2, c2);
waveEpochMovie(c2, F.gm, F.t, 2, files{2}, 'z_band',[5 40], ...
    'Style','phase+event', 'Events',F.E, 'WinMs',[-100 300], 'FPS',o.FPS, ...
    'Trigger','flash', 'EpochLabel','Supplementary Video');
end
