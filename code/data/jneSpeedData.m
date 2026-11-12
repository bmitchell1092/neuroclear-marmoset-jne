function S = jneSpeedData(varargin)
% jneSpeedData  Fig 8c: the speed of every detected wave, spontaneous vs flash-evoked.
%
% Each event is detected with waveEvents and its speed measured with waveEventWavelength
% (temporal frequency x wavelength on the detector's own phase field, at the detected frame),
% the same way for both surfaces:
%   spont   data/fig8/lightRest_*.mat (search window = the window minus 50 ms at each end)
%   flash   data/fig8/flash_*.mat (search window -250..350 ms); windows 'pre' = -150..-50 ms,
%           'evoked' = 110..130 ms (the peak of the event rate)
%
% Speed rises with amplitude (more noise in the phase gradient shortens the apparent
% wavelength), so S.match also gives the spontaneous speeds reweighted to each flash
% window's amplitude deciles.
%
% Output S.E (surface, window, session, speed_m_s, amp, ...)  S.match  S.T (summary)
q = inputParser;
q.addParameter('RhoThr', 0.37);
q.parse(varargin{:}); o = q.Results;
persistent memo memoKey                                % computed once per MATLAB session
if isequal(memoKey, o), S = memo; return; end
dd = fullfile(jneDataDir(), 'fig8');
E = [oneSurface(dir(fullfile(dd,'lightRest_2*.mat')), "spont", o); ...
     oneSurface(dir(fullfile(dd,'flash_*.mat')), "flash", o)];
E.window = repmat("other", height(E), 1);
E.window(E.surface == "spont") = "spont";
E.window(E.surface == "flash" & E.time_ms >= -150 & E.time_ms <= -50) = "pre";
E.window(E.surface == "flash" & E.time_ms >= 110 & E.time_ms <= 130) = "evoked";
E = E(isfinite(E.speed_m_s) & E.speed_m_s > 0 & E.window ~= "other", :);

sp = E.surface == "spont"; match = struct(); T = table();
for g = ["spont_spont","flash_pre","flash_evoked"]
    pr = split(g, "_"); m = E.surface == pr(1) & E.window == pr(2); v = E.speed_m_s(m);
    row = table(pr(1), pr(2), nnz(m), numel(unique(E.session(m))), median(v), prctile(v,25), ...
        prctile(v,75), mean(v >= 0.1 & v <= 0.8), NaN, NaN, NaN, 'VariableNames', ...
        {'surface','window','nEvents','nSess','median','q25','q75','fracReportedRange', ...
         'spontMatchedMedian','spontMatchedQ25','spontMatchedQ75'});
    if pr(1) == "flash"
        w = decileWeights(E.amp(sp), E.amp(m));
        match.(g) = struct('w',w, 'med',wq(E.speed_m_s(sp), w, 0.5), ...
            'iqr',[wq(E.speed_m_s(sp), w, 0.25) wq(E.speed_m_s(sp), w, 0.75)]);
        row.spontMatchedMedian = match.(g).med; row.spontMatchedQ25 = match.(g).iqr(1);
        row.spontMatchedQ75 = match.(g).iqr(2);
    end
    T = [T; row]; %#ok<AGROW>
end
S = struct('E',E, 'match',match, 'T',T, 'params',o);
memo = S; memoKey = o;
end

% =====================================================================
function E = oneSurface(files, surf, o)
[~, ix] = sort({files.name}); files = files(ix);
E = table();
for k = 1:numel(files)
    L = load(fullfile(files(k).folder, files(k).name), 'cube','gridMask','time_ms','session','pitch_mm');
    t = L.time_ms(:).';
    if surf == "spont", dw = [min(t)+50, max(t)-50]; else, dw = [-250 350]; end
    E0 = waveEvents(L.cube, L.gridMask, t, 'RhoThr',o.RhoThr, 'detect_win_ms',dw, 'Verbose',false);
    E1 = waveEventWavelength(E0(:, {'epoch','time_ms','rho','srcRow','srcCol'}), L.cube, L.gridMask, t, ...
            'PixelSpacingMm',L.pitch_mm, 'Verbose',false);
    n = height(E1);
    E = [E; table(repmat(surf,n,1), repmat(string(L.session),n,1), E1.epoch, E1.time_ms, E1.rho, ...
        E1.speed_mm_s/1000, E1.lambda_mm, E1.freqHz, E1.amp, 'VariableNames', {'surface','session', ...
        'epoch','time_ms','rho','speed_m_s','lambda_mm','freqHz','amp'})]; %#ok<AGROW>
    fprintf('[jneSpeedData] %-5s %-34s %5d events\n', surf, L.session, n);
end
end

function w = decileWeights(aS, aT)
% spontaneous events weighted so their amplitude distribution matches the target's
ed = prctile(aT, 0:10:100); ed([1 end]) = [-inf inf]; w = zeros(size(aS));
for b = 1:10
    in = aS >= ed(b) & aS < ed(b+1); if ~any(in), continue; end
    w(in) = mean(aT >= ed(b) & aT < ed(b+1)) / nnz(in);
end
w = w / sum(w);
end

function v = wq(x, w, p)
[xs, ix] = sort(x); cw = cumsum(w(ix)); v = xs(find(cw >= p, 1));
end
