function F = jneFlashWaveData(varargin)
% jneFlashWaveData  Fig 8b: detected wave events per trial around the flash, and the ERP.
%
% Every data/fig8/flash_*.mat cube goes through waveEvents exactly as flashWaveEvents_session
% runs it (rho >= 0.37, search window -250..350 ms). The rate is events per trial per bin,
% averaged over sessions; significance is cluster-corrected against a circular-shift null with
% the session as the unit (waveTimecoursePermTest).
%
% Name/value: 'Win' [-200 350]  'BinMs' 10  'NPerm' 2000  'RhoThr' 0.37  'DetectWinMs' [-250 350]
% Output F.rate (.t .perSession .mean .sem .nSess .nTrials .clusters .sessions)  F.erp (.t .perSession)
q = inputParser;
q.addParameter('Win', [-200 350]);
q.addParameter('BinMs', 10);
q.addParameter('NPerm', 2000);
q.addParameter('RhoThr', 0.37);
q.addParameter('DetectWinMs', [-250 350]);
q.parse(varargin{:}); o = q.Results;
persistent memo memoKey                                % computed once per MATLAB session
if isequal(memoKey, o), F = memo; return; end

files = dir(fullfile(jneDataDir(), 'fig8', 'flash_*.mat'));
[~, ix] = sort({files.name}); files = files(ix);               % by date
edges = o.Win(1):o.BinMs:o.Win(2); tc = edges(1:end-1) + o.BinMs/2;
Y = []; ERP = []; nTr = []; names = strings(0); EV = table(); EP = table();
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'cube','gridMask','time_ms','session');
    cube = S.cube; t = S.time_ms(:).'; n = size(cube,1); nm = string(S.session);
    E = waveEvents(cube, S.gridMask, t, 'RhoThr',o.RhoThr, 'detect_win_ms',o.DetectWinMs, 'Verbose',false);
    Y(end+1,:) = histcounts(E.time_ms, edges) / n; %#ok<AGROW>
    ERP = cat(4, ERP, single(squeeze(mean(double(cube), 1, 'omitnan'))));
    nTr(end+1) = n; names(end+1) = nm; %#ok<AGROW>
    EV = [EV; table(repmat(nm,height(E),1), E.epoch, E.time_ms, 'VariableNames',{'session','epoch','time_ms'})]; %#ok<AGROW>
    EP = [EP; table(repmat(nm,n,1), (1:n)', 'VariableNames',{'session','epoch'})]; %#ok<AGROW>
    fprintf('[jneFlashWaveData] %-34s %4d trials %5d events\n', nm, n, height(E));
end
T = waveTimecoursePermTest(EV, EP, 'Stat','rate', 'Win',o.Win, 'BinMs',o.BinMs, ...
        'NPerm',o.NPerm, 'Verbose',false);
rate = struct('t',tc, 'perSession',Y, 'mean',mean(Y,1), 'sem',std(Y,0,1)/sqrt(size(Y,1)), ...
    'nSess',numel(names), 'nTrials',sum(nTr), 'clusters',T.clusters, 'sessions',names);
F = struct('rate',rate, 'erp',struct('t',t, 'perSession',ERP), 'params',o);
memo = F; memoKey = o;
end
