function C = flashERP_classify(respTr, baseTr, ampP2P, opt)
% flashERP_classify
%
% PURPOSE
% -------
% Single source of truth for the two flash-ERP contact metrics. Given the
% per-flash response- and baseline-window snippets, it computes BOTH axes that
% the 251211 investigation showed are independent:
%
%   1) RESPONSIVENESS (electrode health) - is there a flash-evoked response on
%      this contact at all? Measured on the TRIAL-AVERAGED ERP via a bootstrap
%      CI, so it is robust to single-trial noise: a contact is responsive if the
%      95% CI of (averaged response-window p2p - averaged baseline-window p2p),
%      resampled across flashes, lies above 0, AND the averaged evoked p2p
%      clears AmpFloorUV. This is what tracks contact loss over sessions.
%
%   2) SINGLE-TRIAL DETECTABILITY (SNR) - how reliably does the evoked response
%      beat ongoing LFP on a per-flash basis? This is the old paired signed-rank
%      criterion (now reported as its own axis), plus two continuous summaries:
%        stHitRate : fraction of flashes whose response-window p2p exceeds the
%                    same flash's pre-onset p2p (0.5 = chance)
%        stSNR     : median single-flash response/baseline p2p ratio
%      251211 collapsed on THIS axis (hit-rate ~ chance) while staying healthy
%      on axis 1.
%
% INPUTS
% ------
% respTr : nChan x nRespSamp x nTr  baseline-corrected response-window snippets
% baseTr : nChan x nBaseSamp x nTr  equal-length pre-onset snippets
% ampP2P : nChan x 1                 p2p of the trial-averaged ERP in the
%                                    response window (used for the amplitude floor)
% opt    : struct with fields
%            AmpFloorUV       amplitude floor for both axes (uV)
%            NBoot            bootstrap resamples for responsiveness CI
%            Alpha            signed-rank significance for detectability
%            MinTrialsForTest min flashes to attempt the per-trial test
%            BootSeed         rng seed for the bootstrap ([] = don't reseed)
%
% OUTPUT (C struct, all nChan x 1 unless noted)
% ---------------------------------------------
% C.responsive   logical  robust (averaged-ERP) responsiveness
% C.respCIlo     double   bootstrap 2.5% CI of (respP2P-baseP2P) on the mean ERP
% C.baseP2P_avg  double   p2p of the trial-averaged ERP in the baseline window
% C.detectable   logical  single-trial signed-rank detectability
% C.pval         double   paired signed-rank p (per-trial resp vs base p2p)
% C.stHitRate    double   fraction of flashes with resp p2p > base p2p
% C.stSNR        double   median single-flash resp/base p2p ratio
%
% Blake / U01 pipeline

nChan = size(respTr, 1);
nTr   = size(respTr, 3);

% Per-trial p2p in each window (nChan x nTr)
respP2P_tr = squeeze(max(respTr, [], 2) - min(respTr, [], 2));
baseP2P_tr = squeeze(max(baseTr, [], 2) - min(baseTr, [], 2));
if nTr == 1, respP2P_tr = respP2P_tr(:); baseP2P_tr = baseP2P_tr(:); end

% Averaged-ERP p2p in each window
p2pAvg = @(M) max(M, [], 2) - min(M, [], 2);
baseP2P_avg = p2pAvg(mean(baseTr, 3));

% ---------------- Axis 1: robust responsiveness (bootstrap on averaged ERP) ----
respCIlo = nan(nChan, 1);          % stays NaN in the degenerate (1-flash) case
floorOK  = ampP2P >= opt.AmpFloorUV;

if nTr >= 2
    if ~isempty(opt.BootSeed), rng(opt.BootSeed); end
    B = opt.NBoot;
    diffsBoot = zeros(nChan, B);
    for b = 1:B
        idx = randi(nTr, nTr, 1);
        mr  = mean(respTr(:, :, idx), 3);     % nChan x nRespSamp
        mb  = mean(baseTr(:, :, idx), 3);
        diffsBoot(:, b) = p2pAvg(mr) - p2pAvg(mb);
    end
    respCIlo   = quantile(diffsBoot, 0.025, 2);
    responsive = (respCIlo > 0) & floorOK;
else
    % Degenerate (single flash): fall back to a point estimate.
    responsive = (p2pAvg(mean(respTr,3)) > baseP2P_avg) & floorOK;
end

% ---------------- Axis 2: single-trial detectability ----------------------------
stHitRate = mean(respP2P_tr > baseP2P_tr, 2);
stSNR     = median(respP2P_tr ./ (baseP2P_tr + eps), 2);

pval       = nan(nChan, 1);
detectable = false(nChan, 1);
if nTr >= opt.MinTrialsForTest
    for c = 1:nChan
        d = respP2P_tr(c, :) - baseP2P_tr(c, :);
        if any(d ~= 0)
            pval(c) = signrank(respP2P_tr(c, :), baseP2P_tr(c, :));
        else
            pval(c) = 1;
        end
        detectable(c) = isfinite(pval(c)) && pval(c) < opt.Alpha ...
            && median(d) > 0 && floorOK(c);
    end
end

% ---------------- Package ----------------
C = struct();
C.responsive  = responsive;
C.respCIlo    = respCIlo;
C.baseP2P_avg = baseP2P_avg;
C.detectable  = detectable;
C.pval        = pval;
C.stHitRate   = stHitRate;
C.stSNR       = stSNR;
end
