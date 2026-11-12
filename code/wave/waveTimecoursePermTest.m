function T = waveTimecoursePermTest(EV, EP, varargin)
% waveTimecoursePermTest   Which parts of a wave-event time course are real?
%                          Surrogate triggers + cluster-based permutation.
%
% PURPOSE
% -------
% A binned rate or alignment curve wiggles. Reading a peak off one is how a window
% gets chosen, and once a window is chosen from the data the test at that window is
% circular. This replaces window-picking with a statement of the form "the curve
% departs from the trigger-free null over THIS contiguous interval, corrected for
% having looked at every bin".
%
% THE NULL: CIRCULAR SHIFT WITHIN EPOCH
% -------------------------------------
% Each epoch's event train is rolled by its own random offset inside the detection
% window. That keeps everything except the thing being tested:
%   * each epoch keeps its exact event COUNT (so the rate null has the right spread
%     across sessions, which differ several-fold in baseline rate);
%   * each epoch keeps its own inter-event structure;
%   * the pooled DIRECTION distribution is untouched, so an R null reflects the
%     alignment obtainable by chance from this array's standing direction bias;
%   * only the relationship to t=0 is destroyed.
% A flat rate null is therefore the correct expectation, and any structure in the
% observed curve is trigger-locked by construction.
%
% This is the cheap surrogate. The stronger one -- re-drawing TRIGGER TIMES at random
% and rebuilding cubes -- also resamples which LFP the detector saw, and would catch
% an artifact where the epoch selection itself (isolation, eyes-open gating) imposes
% temporal structure. It costs a full recompute per permutation and is not done here.
% Treat a cluster that survives this test as real time-locking, not yet as proof that
% the epoch selection is innocent.
%
% SESSION IS THE UNIT. The statistic is the across-session mean of a per-session
% curve, and every permutation re-derives it the same way. Pooling events instead
% would let the largest sessions carry the result.
%
% CLUSTERING. Per bin, the observed value is z-scored against that bin's own null
% distribution (so a bin with fewer events is not automatically more extreme).
% Contiguous runs of |z| > ClusterZ form clusters; a cluster's mass is the sum of its
% z. Each permutation contributes its own largest cluster mass, and a cluster is
% significant when it exceeds the (1-Alpha) quantile of that max-cluster
% distribution. That correction is what makes "there is a bump at 115 ms" a claim
% about the whole curve rather than about the bin someone looked at.
%
% INPUTS
% ------
% EV  events table: .session, .epoch, the time column named by 'TimeCol', and
%     .dirDeg when Stat is 'R'
% EP  epoch table:  .session, .epoch  (one row per epoch -- the rate denominator)
% Name/value:
%   'Stat'       'rate' | 'R' | 'polarity' | 'frac'   quantity to test
%   'AttrCol'    ''     REQUIRED for 'frac': name of a NUMERIC per-event column of EV.
%                       The curve is that column's mean over the events in each bin
%                       (per session, then across sessions) -- so it is a quantity
%                       CONDITIONAL ON AN EVENT HAVING HAPPENED, and is by construction
%                       independent of the event rate.
%                       ⭐ This is what the motif analysis needs and why it exists. A
%                       COUNT of motif events per bin rises whenever the rate rises, and
%                       the saccade rate rise is already established -- so a count can
%                       never distinguish "the saccade made a pattern" from "the saccade
%                       made more waves". With AttrCol = a 0/1 motif-membership flag, the
%                       curve is the FRACTION of waves that are stereotyped, which can
%                       fall while the rate rises. Use it with a continuous stereotypy
%                       score too (mean pattern similarity), which needs no clustering
%                       threshold at all.
%   'RefDirDeg'  NaN    REQUIRED for 'polarity': the reference sense (top-CW deg).
%                       polarity = kernel-weighted mean of cos(theta - ref), in [-1 1].
%                       +1 = every wave travels toward ref, 0 = symmetric, -1 = away.
%                       ⚠️ AXIS statistics are mod-180 and BLIND to this by construction:
%                       a 3:1 polarised distribution and a symmetric one give identical
%                       R values. Polarity is the mod-360 question and needs its own stat.
%   'TimeCol'    'time_ms'      which time base (e.g. 't_off' for landing-locked)
%   'Win'        [-250 350]     detection window; the shift is circular INSIDE this
%   'BinMs'      10
%   'NPerm'      2000
%   'ClusterZ'   2
%   'Alpha'      0.05
%   'MinEvBin'   3              bins with fewer events give no R (left NaN)
%   'Seed'       7
%   'Verbose'    true
%
% OUTPUT (struct T)
%   .t .obs .nullMean .nullSd .z          per bin
%   .nullMat  NPerm x nBin surrogate curves -- kept so a caller can build a null for
%             a DERIVED quantity (e.g. the spectrum of the curve) instead of
%             inventing a second, differently-constructed surrogate for it
%   .clusters  table: t1, t2, mass, p, sign
%   .Params
%
% EXAMPLE
%   T = waveTimecoursePermTest(EV, EP, 'Stat','R', 'TimeCol','t_off');
%   T.clusters(T.clusters.p < 0.05, :)
%
% See also: waveEvents, run_saccadeWaveTimecourse, waveEventRateTest

p = struct('Stat','rate', 'RefDirDeg',NaN, 'AttrCol','', 'TimeCol','time_ms', 'Win',[-250 350], 'BinMs',10, ...
           'Estimator','bin', 'BandwidthMs',[], 'EvalStepMs',[], 'FineMs',1, ...
           'EdgeSigma',3, 'MinNeff',3, ...
           'NPerm',2000, 'ClusterZ',2, 'Alpha',0.05, 'MinEvBin',3, 'Seed',7, ...
           'Verbose',true);
for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end
rng(p.Seed);
useKernel = strcmpi(p.Estimator,'kernel');
if useKernel && isempty(p.BandwidthMs)
    error('waveTimecoursePermTest:noBandwidth', ...
        'Estimator ''kernel'' requires BandwidthMs. Choose it from n_eff, not by eye.');
end
if isempty(p.EvalStepMs), p.EvalStepMs = p.BinMs; end

span = diff(p.Win);
if ~useKernel
    edges = p.Win(1):p.BinMs:p.Win(2);
    ctr   = edges(1:end-1) + p.BinMs/2;
    K = []; K2 = []; evalIdx = [];
else
    % ---- kernel evaluation grid ----
    % EDGE TRIM. A Gaussian centred within EdgeSigma sigma of the window edge collects
    % weight from outside the data range, so its density is biased DOWN -- the same
    % class of error as the landing-locked coverage taper, and it would read as a
    % suppression at both ends. The curve is therefore only evaluated where the kernel
    % is fully supported.
    pad  = p.EdgeSigma * p.BandwidthMs;
    ctr  = (p.Win(1)+pad) : p.EvalStepMs : (p.Win(2)-pad);
    if numel(ctr) < 5
        error('waveTimecoursePermTest:bandwidthTooWide', ...
            ['BandwidthMs %g leaves only %d evaluation points inside Win [%g %g] ' ...
             'after the %g-sigma edge trim.'], p.BandwidthMs, numel(ctr), p.Win, p.EdgeSigma);
    end
    % ---- fine lattice + kernel, so the estimate is a CONVOLUTION ----
    % Evaluating exp() per event per eval point would be ~1e9 operations across 2000
    % permutations. Binning finely (FineMs, negligible against sigma) and convolving is
    % mathematically the same kernel estimate and runs in milliseconds. sumW comes from
    % convolving counts with K; sumW2 from convolving with K.^2 -- which is what the
    % effective sample size needs and is easy to get wrong by reusing K twice.
    fineEdges = p.Win(1):p.FineMs:p.Win(2);
    fineCtr   = fineEdges(1:end-1) + p.FineMs/2;
    kt = -ceil(4*p.BandwidthMs/p.FineMs)*p.FineMs : p.FineMs : ...
          ceil(4*p.BandwidthMs/p.FineMs)*p.FineMs;
    K  = exp(-kt.^2 / (2*p.BandwidthMs^2)) / (p.BandwidthMs*sqrt(2*pi));
    K  = K / (sum(K)*p.FineMs);         % normalise to a density over ms
    K2 = K.^2;
    edges = fineEdges;
    [~, evalIdx] = min(abs(fineCtr(:) - ctr(:).'), [], 1);
end
nB = numel(ctr);

% ---- index events and epochs by session ----
us  = unique(EP.session, 'stable');
nS  = numel(us);
[~, sEv] = ismember(EV.session, us);
[~, sEp] = ismember(EP.session, us);
keep = sEv > 0;
t0   = EV.(p.TimeCol)(keep);
sEv  = sEv(keep);
% dEv is the generic PER-EVENT PAYLOAD carried through the shift with its event:
% direction in radians for 'R'/'polarity', the AttrCol value for 'frac', unused for 'rate'.
dEv  = zeros(size(t0));
if any(strcmpi(p.Stat,{'R','polarity'})), dEv = deg2rad(EV.dirDeg(keep)); end
if strcmpi(p.Stat,'frac')
    if isempty(p.AttrCol) || ~ismember(p.AttrCol, EV.Properties.VariableNames)
        error('waveTimecoursePermTest:noAttr', ...
            "Stat 'frac' needs 'AttrCol' naming a numeric column of EV (got '%s').", p.AttrCol);
    end
    dEv = double(EV.(p.AttrCol)(keep));
    if all(isnan(dEv))
        error('waveTimecoursePermTest:attrAllNaN','EV.%s is entirely NaN.', p.AttrCol);
    end
end
if strcmpi(p.Stat,'polarity')
    assert(isfinite(p.RefDirDeg), 'waveTimecoursePermTest:noRef', ...
        "Stat 'polarity' needs 'RefDirDeg' -- polarity is meaningless without a reference sense.");
end

% Epoch identity within session, so the roll is per epoch and not per session.
epKey = EV.session(keep) + "|" + string(EV.epoch(keep));
[~, ~, epIdx] = unique(epKey);
nEpAll = max(epIdx);

inWin = t0 >= p.Win(1) & t0 <= p.Win(2);
t0 = t0(inWin); sEv = sEv(inWin); dEv = dEv(inWin); epIdx = epIdx(inWin);

nEpSess = accumarray(sEp, 1, [nS 1]);
if p.Verbose
    fprintf('waveTimecoursePermTest[%s|%s]: %d events, %d epochs, %d sessions, %d perms\n', ...
        p.Stat, p.TimeCol, numel(t0), sum(nEpSess), nS, p.NPerm);
end

% ---- observed ----
obs = statCurve(t0, sEv, dEv, edges, nS, nB, nEpSess, p, K, K2, evalIdx);

% ---- surrogates ----
nullMat = nan(p.NPerm, nB);
for k = 1:p.NPerm
    shift = rand(nEpAll,1) * span;                 % one offset per EPOCH
    ts = mod(t0 - p.Win(1) + shift(epIdx), span) + p.Win(1);
    nullMat(k,:) = statCurve(ts, sEv, dEv, edges, nS, nB, nEpSess, p, K, K2, evalIdx);
end

nullMean = mean(nullMat, 1, 'omitnan');
nullSd   = std(nullMat, 0, 1, 'omitnan');
z        = (obs - nullMean) ./ max(nullSd, eps);

% ---- cluster mass, and the null distribution of the LARGEST cluster ----
obsCl = findClusters(z, p.ClusterZ, ctr);
% Fast path for the null loop: only the LARGEST cluster mass is needed, so it is
% computed numerically instead of building a table 2000 times over.
maxNull = zeros(p.NPerm,1);
for k = 1:p.NPerm
    zk = (nullMat(k,:) - nullMean) ./ max(nullSd, eps);
    maxNull(k) = maxClusterMass(zk, p.ClusterZ);
end
if ~isempty(obsCl)
    obsCl.p = arrayfun(@(m) (1 + sum(maxNull >= abs(m))) / (p.NPerm + 1), obsCl.mass);
    obsCl = sortrows(obsCl, 'p');
end

T = struct('t',ctr, 'obs',obs, 'nullMean',nullMean, 'nullSd',nullSd, 'z',z, ...
           'nullMat',nullMat, ...
           'clusters',obsCl, 'maxNull',maxNull, 'Params',p);

if p.Verbose
    if isempty(obsCl)
        fprintf('  no cluster reached |z| > %g\n', p.ClusterZ);
    else
        sig = obsCl(obsCl.p < p.Alpha, :);
        fprintf('  %d cluster(s), %d significant at alpha=%g\n', height(obsCl), height(sig), p.Alpha);
        for i = 1:height(obsCl)
            fprintf('    %6.0f to %6.0f ms  mass %+8.1f  p = %.4f%s\n', ...
                obsCl.t1(i), obsCl.t2(i), obsCl.mass(i), obsCl.p(i), ...
                repmat(' *', 1, obsCl.p(i) < p.Alpha));
        end
    end
end
end

% =====================================================================
function y = kernelCurve(t, sEv, dEv, edges, nS, ~, nEpSess, p, K, K2, evalIdx)
% Gaussian-kernel estimate, per session, then the across-session mean.
%
% WHY THIS EXISTS. Binning has two different costs and only one of them is obvious.
% For RATE it is just a smoothing choice -- unbiased at any width, trading variance
% against resolution. For ALIGNMENT it is a COMPUTABILITY problem: the resultant needs
% a minimum count per session-bin cell or it is undefined, and it is biased upward at
% small n, so the bin width changes the estimator's BIAS. On the Atlas array only 25%
% of session x 10 ms cells reached that floor, so three-quarters of the curve was NaN.
%
% A kernel dissolves that. Every evaluation point uses ALL of a session's events,
% weighted by their distance in time, so there is no cell membership and no hard
% floor -- the effective sample size becomes continuous:
%
%     n_eff = (sum w)^2 / sum w^2
%
% and the same small-n bias correction applies with n_eff in place of a count. It also
% removes fixed bin EDGES, which were demonstrably moving cluster boundaries by a bin
% under small re-alignments of the trigger.
nF = numel(edges) - 1;
b  = discretize(t, edges);
ok = ~isnan(b);
b = b(ok); s = sEv(ok); d = dEv(ok);
idx = sub2ind([nS nF], s, b);

cnt = reshape(accumarray(idx, 1, [nS*nF 1]), nS, nF);
sw  = conv2(cnt, K,  'same');           % sum of weights            -> density
sw2 = conv2(cnt, K2, 'same');           % sum of SQUARED weights    -> n_eff

if strcmpi(p.Stat,'rate')
    % K is normalised as a density over ms, so sw is events per epoch per ms.
    M = (sw(:, evalIdx) ./ nEpSess) * 1000;
    y = mean(M, 1, 'omitnan');
    return;
end

if strcmpi(p.Stat,'frac')
    % Kernel-weighted mean of the per-event attribute: convolve the attribute SUM with
    % the same K used for the count, then divide. Weighting numerator and denominator
    % identically is what keeps this a mean and not a density.
    dz = d; dz(isnan(dz)) = 0;
    sv = reshape(accumarray(idx, dz, [nS*nF 1]), nS, nF);
    wv = conv2(sv, K, 'same');
    M  = wv(:, evalIdx) ./ max(sw(:, evalIdx), eps);
    M(((sw(:, evalIdx).^2) ./ max(sw2(:, evalIdx), eps)) < p.MinNeff) = NaN;
    y = mean(M, 1, 'omitnan');
    return;
end

cs = reshape(accumarray(idx, cos(d), [nS*nF 1]), nS, nF);
sn = reshape(accumarray(idx, sin(d), [nS*nF 1]), nS, nF);
wc = conv2(cs, K, 'same');
ws = conv2(sn, K, 'same');

if strcmpi(p.Stat,'polarity')
    % Projection of each event's direction onto the reference sense, kernel-weighted.
    % Reuses the same wc/ws convolutions the R path uses, so the two statistics are
    % measured on identical smoothing and the permutation null transfers unchanged.
    rr = deg2rad(p.RefDirDeg);
    M  = (wc(:, evalIdx)*cos(rr) + ws(:, evalIdx)*sin(rr)) ./ max(sw(:, evalIdx), eps);
    M(((sw(:, evalIdx).^2) ./ max(sw2(:, evalIdx), eps)) < p.MinNeff) = NaN;
    y = mean(M, 1, 'omitnan');
    return;
end

swE  = sw(:, evalIdx);  sw2E = sw2(:, evalIdx);
R    = hypot(wc(:, evalIdx), ws(:, evalIdx)) ./ max(swE, eps);
nEff = (swE.^2) ./ max(sw2E, eps);

% Bias correction with n_eff, exactly as the binned path does with a count. Without
% it a sparse session reads as more aligned purely because it has fewer events near
% that time point -- which is the whole failure mode this is meant to avoid.
M = sqrt(max(0, (nEff.*R.^2 - 1) ./ max(nEff - 1, eps)));
M(nEff < p.MinNeff) = NaN;
y = mean(M, 1, 'omitnan');
end

% =====================================================================
function y = statCurve(t, sEv, dEv, edges, nS, nB, nEpSess, p, K, K2, evalIdx)
if strcmpi(p.Estimator,'kernel')
    y = kernelCurve(t, sEv, dEv, edges, nS, nB, nEpSess, p, K, K2, evalIdx);
    return;
end
% Per-session curve, then the across-session mean. Sessions with no events in a bin
% contribute a rate of 0 (real) but no R (undefined) -- the two must not be conflated.
b = discretize(t, edges);
ok = ~isnan(b);
b = b(ok); s = sEv(ok); d = dEv(ok);
idx = sub2ind([nS nB], s, b);

cnt = accumarray(idx, 1, [nS*nB 1]);
cnt = reshape(cnt, nS, nB);

if strcmpi(p.Stat,'rate')
    M = cnt ./ (nEpSess * p.BinMs/1000);
    y = mean(M, 1, 'omitnan');
elseif strcmpi(p.Stat,'frac')
    dz = d; dz(isnan(dz)) = 0;
    sv = reshape(accumarray(idx, dz, [nS*nB 1]), nS, nB);
    M = sv ./ max(cnt, 1);
    M(cnt < p.MinEvBin) = NaN;      % a fraction off 1-2 events is 0 or 1 by construction
    y = mean(M, 1, 'omitnan');
elseif strcmpi(p.Stat,'polarity')
    cs = reshape(accumarray(idx, cos(d), [nS*nB 1]), nS, nB);
    sn = reshape(accumarray(idx, sin(d), [nS*nB 1]), nS, nB);
    rr = deg2rad(p.RefDirDeg);
    M = (cs*cos(rr) + sn*sin(rr)) ./ max(cnt, 1);
    M(cnt < p.MinEvBin) = NaN;
    y = mean(M, 1, 'omitnan');
else
    cs = reshape(accumarray(idx, cos(d), [nS*nB 1]), nS, nB);
    sn = reshape(accumarray(idx, sin(d), [nS*nB 1]), nS, nB);
    M = hypot(cs, sn) ./ max(cnt, 1);
    M(cnt < p.MinEvBin) = NaN;      % R off 1-2 events is 1.0 by construction
    y = mean(M, 1, 'omitnan');
end
end

% =====================================================================
function m = maxClusterMass(z, thr)
% Largest |summed z| over sign-consistent runs above threshold. Same definition as
% findClusters, without the bookkeeping -- these two MUST agree or the null is being
% built for a different statistic than the observed one is measured with.
%
% The body now lives in waveClusterMass.m so that callers building a FAMILY-WISE null on
% top of .nullMat (waveMotifOrganization corrects across motifs as well as time) use this
% exact definition rather than a fourth copy of it.
m = waveClusterMass(z, thr);
end

% =====================================================================
function C = findClusters(z, thr, ctr)
% Contiguous runs of |z| > thr, signed, with summed z as mass.
C = table('Size',[0 4], 'VariableTypes',{'double','double','double','double'}, ...
          'VariableNames',{'t1','t2','mass','sign'});
zz = z; zz(~isfinite(zz)) = 0;
above = abs(zz) > thr;
if ~any(above), return; end
d = diff([false above false]);
i1 = find(d == 1); i2 = find(d == -1) - 1;
for k = 1:numel(i1)
    seg = zz(i1(k):i2(k));
    % A run must not mix signs -- a +z bin adjacent to a -z bin is two effects, and
    % summing them would cancel a real suppression against a real enhancement.
    if all(seg > 0) || all(seg < 0)
        C = [C; {ctr(i1(k)), ctr(i2(k)), sum(seg), sign(seg(1))}]; %#ok<AGROW>
    else
        pos = seg > 0;
        C = [C; splitRun(seg, pos,  ctr, i1(k))]; %#ok<AGROW>
        C = [C; splitRun(seg, ~pos, ctr, i1(k))]; %#ok<AGROW>
    end
end
end

function C = splitRun(seg, mask, ctr, off)
C = table('Size',[0 4], 'VariableTypes',{'double','double','double','double'}, ...
          'VariableNames',{'t1','t2','mass','sign'});
d = diff([false mask(:).' false]);
i1 = find(d == 1); i2 = find(d == -1) - 1;
for k = 1:numel(i1)
    s = seg(i1(k):i2(k));
    C = [C; {ctr(off+i1(k)-1), ctr(off+i2(k)-1), sum(s), sign(s(1))}]; %#ok<AGROW>
end
end
