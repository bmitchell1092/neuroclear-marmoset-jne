function R = erpFromTrials(Tr, tMs, channelMap, opt)
% erpFromTrials  Flash ERP metrics from the per-trial tensor (Fig 7).
%
% Tr          nChan x nT x nTrials, uV, already baseline-corrected (as exported)
% tMs         1 x nT, ms from flash onset
% channelMap  8 x 10, 0-based channel per grid position (NaN = no contact)
% opt         the flashERP_session options saved with the data (BaselineMs, RespMs,
%             AmpFloorUV, NBoot, Alpha, MinTrialsForTest, BootSeed)
%
% This is the part of flashERP_session (wave-DAQ) that runs after the trials are cut,
% copied line for line: ERP, peak-to-peak in RespMs, peak latency, the responsive /
% detectable classification (flashERP_classify) and the grid maps.
nChan = size(Tr,1); nTr = size(Tr,3);
baseSel     = tMs >= opt.BaselineMs(1) & tMs <= opt.BaselineMs(2);
respSel     = tMs >= opt.RespMs(1)     & tMs <= opt.RespMs(2);
respLen     = opt.RespMs(2) - opt.RespMs(1);
baseTestSel = tMs >= -respLen & tMs < 0;

erp = mean(Tr, 3);
sem = std(Tr, 0, 3) / sqrt(nTr);
p2p = @(v) max(v, [], 2) - min(v, [], 2);
ampP2P      = p2p(erp(:, respSel));
baseP2P_erp = p2p(erp(:, baseTestSel));
snr         = ampP2P ./ (baseP2P_erp + eps);
baselineRMS = sqrt(mean(erp(:, baseSel).^2, 2));

peakLatencyMs = nan(nChan,1);
respMsVals = tMs(respSel);
[~, iPk] = max(abs(erp(:, respSel)), [], 2);
peakLatencyMs(:) = respMsVals(iPk);

C = flashERP_classify(Tr(:, respSel, :), Tr(:, baseTestSel, :), ampP2P, ...
    struct('AmpFloorUV', opt.AmpFloorUV, 'NBoot', opt.NBoot, ...
           'Alpha', opt.Alpha, 'MinTrialsForTest', opt.MinTrialsForTest, ...
           'BootSeed', opt.BootSeed));
responsive = C.responsive; detectable = C.detectable;

[nr, nc] = size(channelMap); chMap1 = channelMap + 1;
f = {'ampP2P','responsive','detectable','stHitRate','stSNR','snr','peakLatencyMs'};
v = {ampP2P, double(responsive), double(detectable), C.stHitRate, C.stSNR, snr, peakLatencyMs};
map = struct();
for i = 1:numel(f)
    m = nan(nr, nc);
    for rr = 1:nr
        for cc = 1:nc
            ch = chMap1(rr,cc);
            if ~isnan(ch) && ch >= 1 && ch <= nChan, m(rr,cc) = v{i}(ch); end
        end
    end
    map.(f{i}) = m;
end

R = struct('nChan',nChan, 'nTrials',nTr, 'tMs',tMs, 'erp',erp, 'sem',sem, 'ampP2P',ampP2P, ...
    'snr',snr, 'pval',C.pval, 'responsive',responsive, 'respCIlo',C.respCIlo, ...
    'detectable',detectable, 'stHitRate',C.stHitRate, 'stSNR',C.stSNR, ...
    'peakLatencyMs',peakLatencyMs, 'baselineRMS',baselineRMS, ...
    'nResponsive',sum(responsive), 'respFrac',mean(responsive), ...
    'channelMap',channelMap, 'map',map, 'params',opt, 'trials',Tr);
end
