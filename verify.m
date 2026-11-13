function ok = verify()
% verify  Recompute the numbers behind Figures 7 and 8 from data/ and check them against
% data/expected.mat (the values in the published figures). Prints one line per check and
% the main numbers quoted in the paper.
setup_paths();
X = load(fullfile(jneDataDir(), 'expected.mat'));
res = {};
tol = 1e-9;     % Fig 7, and anything that is a count
tolW = 1e-6;    % wave rho and speed: small floating-point differences across MATLAB versions or
                % machines are tolerated; counts of trials and events must match exactly
chk = @(name, pass, detail) fprintf('  %-4s %-44s %s\n', ternary(pass,'ok','FAIL'), name, detail);

fprintf('\nFigure 7\n');
D = jneData('Verbose', false);
Qx = X.fig7.Q; [~, ia] = ismember(Qx.Date, D.Q.Date);
p = all(ia > 0) && isequal(D.Q.nResponsive(ia), Qx.nResponsive);
chk('responsive contacts per session (g)', p, mat2str(D.Q.nResponsive(ia).')); res{end+1} = p;
for i = 1:numel(X.fig7.rep)
    e = X.fig7.rep(i); r = D.flash(string({D.flash.date}) == e.date).R;
    d = max(abs(r.erp - e.erp), [], 'all');
    p = d < tol && isequal(r.responsive(:), e.responsive(:)) && max(abs(r.ampP2P - e.ampP2P)) < tol ...
        && isequal(r.peakLatencyMs, e.peakLatencyMs) && max(abs(r.sem - e.sem), [], 'all') < tol;
    chk(sprintf('ERP, p2p, responsive, latency, SEM (%s)', e.date), p, sprintf('max ERP diff %.1e uV', d));
    res{end+1} = p; %#ok<AGROW>
end
p = all(arrayfun(@(k) isequal(D.darkPSD(k).Pxx, X.fig7.psd(k).Pxx), 1:numel(D.darkPSD)));
chk('resting spectra (f)', p, ''); res{end+1} = p;

fprintf('\nFigure 8\n');
F = jneFlashWaveData(); x = X.fig8.rate;
p = F.rate.nSess == x.nSess && F.rate.nTrials == x.nTrials;
chk('flash sessions / trials (b)', p, sprintf('%d / %d', F.rate.nSess, F.rate.nTrials)); res{end+1} = p;
d = max(abs(F.rate.perSession - x.perSession), [], 'all');
chk('wave events per trial per bin (b)', d < tol, sprintf('max diff %.1e', d)); res{end+1} = d < tol;
c1 = F.rate.clusters; c2 = x.clusters;
p = height(c1) == height(c2) && isequal(c1{:,{'t1','t2','sign'}}, c2{:,{'t1','t2','sign'}}) && max(abs(c1.p - c2.p)) < tol;
chk('significant clusters (b)', p, strjoin(compose('%d-%d ms', c1.t1(c1.p<0.05 & c1.sign>0), c1.t2(c1.p<0.05 & c1.sign>0)), ', '));
res{end+1} = p;
erp = squeeze(mean(F.erp.perSession, [1 2 4], 'omitnan')); d = max(abs(erp - X.fig8.erpMean));
chk('evoked LFP, array mean (b)', d < 1e-4, sprintf('max diff %.1e uV', d)); res{end+1} = d < 1e-4;
S = jneSpeedData(); E = S.E;
vs = sort(E.speed_m_s(E.surface == "spont")); ve = sort(E.speed_m_s(E.window == "evoked"));
p = isequal(size(vs), size(X.fig8.speedSpont)) && max(abs(vs - X.fig8.speedSpont)) < tolW;
chk('spontaneous wave speeds (c)', p, sprintf('n = %d, median %.2f m/s', numel(vs), median(vs))); res{end+1} = p;
p = isequal(size(ve), size(X.fig8.speedEvoked)) && max(abs(ve - X.fig8.speedEvoked)) < tolW;
chk('flash-evoked wave speeds (c)', p, sprintf('n = %d, median %.2f m/s', numel(ve), median(ve))); res{end+1} = p;
Ps = jneSpontEpochData(); Pf = jneFlashEpochData();
p = sameEvents(Ps.E, X.fig8.spontEvents);
chk('spontaneous exemplar detections (a)', p, sprintf('%d events', height(Ps.E))); res{end+1} = p;
p = sameEvents(Pf.E, X.fig8.flashEvents);
chk('flash exemplar detections (a)', p, sprintf('%d events', height(Pf.E))); res{end+1} = p;

ok = all([res{:}]);
fprintf('\n%d of %d checks passed.\n', nnz([res{:}]), numel(res));

fprintf('\nNumbers in the text and legends\n');
fprintf('  Fig 8b: %d sessions, %d trials\n', F.rate.nSess, F.rate.nTrials);
T = S.T;
for k = 1:height(T)
    fprintf('  Fig 8c %-6s %-7s n = %5d, median %.2f m/s (IQR %.2f-%.2f), %.0f%% in 0.1-0.8 m/s', ...
        T.surface(k), T.window(k), T.nEvents(k), T.median(k), T.q25(k), T.q75(k), 100*T.fracReportedRange(k));
    if isfinite(T.spontMatchedMedian(k))
        fprintf('; spontaneous at this amplitude %.2f m/s', T.spontMatchedMedian(k));
    end
    fprintf('\n');
end
end

function p = sameEvents(a, b)
v = {'epoch','time_ms','rho','srcRow','srcCol'};
p = height(a) == height(b) && max(abs(a{:,v} - b{:,v}), [], 'all') < 1e-6;
if height(a) == 0 && height(b) == 0, p = true; end
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
