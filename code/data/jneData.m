function D = jneData(varargin)
% jneData  Everything Figure 7 draws, built from data/fig7.
%
%   D.flash    the three representative sessions (2 wk / 2 mo / 4 mo): erpFromTrials result,
%              slimmed, plus the QC good-channel mask
%   D.other    the other Fig 7 sessions (1 wk, ~7 wk), same form (panel c uses 1 wk)
%   D.Q        responsive-contact yield per session (panel g)
%   D.darkPSD  resting spectra for the three representative dates (panel f)
%
% Name/value: 'Sessions' ["251028" "251219" "260203"]  'ImplantDate' '2025-10-14'
p = inputParser;
p.addParameter('Sessions', ["251028" "251219" "260203"]);
p.addParameter('ImplantDate', '2025-10-14');
p.addParameter('Verbose', true);
p.parse(varargin{:}); opt = p.Results;
implant = datetime(opt.ImplantDate);
rep = string(opt.Sessions);
dd = fullfile(jneDataDir(), 'fig7');

files = dir(fullfile(dd, 'flash_*.mat'));
Q = table(); flash = struct('date',{},'ageDays',{},'label',{},'R',{}); other = flash;
for i = 1:numel(files)
    S = load(fullfile(dd, files(i).name));
    if opt.Verbose, fprintf('[jneData] %s: %d trials\n', S.date, size(S.trials,3)); end
    R = erpFromTrials(S.trials, S.tMs, S.channelMap, S.params);
    R.expID = S.session;
    age = days(datetime(char(S.date), 'InputFormat','yyMMdd') - implant);
    Q = [Q; table(S.date, age, R.nResponsive, R.respFrac, R.nChan, ...
        'VariableNames', {'Date','ageDays','nResponsive','respFrac','nChan'})]; %#ok<AGROW>
    slim = slimR(R); slim.goodMask = S.goodMask(:);
    e = struct('date',S.date, 'ageDays',age, 'label',ageLabel(age), 'R',slim);
    if any(S.date == rep), flash(end+1) = e; else, other(end+1) = e; end %#ok<AGROW>
end
[~, o] = sort(Q.ageDays); Q = Q(o,:);
[~, loc] = ismember(rep, string({flash.date})); flash = flash(loc);

darkPSD = struct('date',{},'ageDays',{},'label',{},'f',{},'Pxx',{});
for k = 1:numel(rep)
    S = load(fullfile(dd, sprintf('darkRest_psd_%s.mat', rep(k))));
    age = days(datetime(char(rep(k)), 'InputFormat','yyMMdd') - implant);
    darkPSD(end+1) = struct('date',rep(k), 'ageDays',age, 'label',ageLabel(age), ...
        'f',S.f, 'Pxx',S.Pxx); %#ok<AGROW>
end

D = struct('flash',flash, 'other',other, 'Q',Q, 'darkPSD',darkPSD, 'implant',implant, ...
    'sessions',rep, 'channelMap',flash(1).R.channelMap, 'builtFrom',string(opt.ImplantDate));
end

% =====================================================================
function r = slimR(R)
r = struct('expID',R.expID, 'tMs',R.tMs, 'erp',R.erp, 'ampP2P',R.ampP2P, ...
    'peakLatencyMs',R.peakLatencyMs, 'snr',R.snr, 'responsive',R.responsive(:), ...
    'channelMap',R.channelMap, 'map',R.map, 'RespMs',R.params.RespMs, 'nChan',R.nChan, ...
    'nTrials',R.nTrials, 'sem',std(R.trials, 0, 3, 'omitnan') ./ sqrt(size(R.trials,3)));
end

function lab = ageLabel(ageDays)
if ageDays < 56, lab = sprintf('%d wk', max(0,round(ageDays/7)));
else, lab = sprintf('%d mo', round(ageDays/30.4)); end
end
