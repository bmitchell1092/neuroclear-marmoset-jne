function m = waveClusterMass(z, thr)
% waveClusterMass   Largest |summed z| over sign-consistent runs above threshold.
%
% PURPOSE
% -------
% The cluster statistic, in ONE place. `waveTimecoursePermTest` builds its null from this
% and measures its observed curve with the same definition; anything that builds a
% FAMILY-WISE null on top (e.g. `waveMotifOrganization`, correcting across motifs as well
% as across time) must use the identical function or the two nulls describe different
% statistics. That warning was already written in waveTimecoursePermTest's local copy --
% "these two MUST agree or the null is being built for a different statistic than the
% observed one is measured with" -- so the copy became a file rather than being duplicated
% a third time.
%
% A run must not mix signs: a +z bin adjacent to a -z bin is two effects, not one big one,
% so each sign is scanned separately.
%
% INPUTS   z    1 x nBin z-scores (NaN treated as 0)
%          thr  cluster-forming threshold on |z|
% OUTPUT   m    largest cluster mass (>= 0); 0 when no bin clears thr
%
% See also: waveTimecoursePermTest, waveMotifOrganization

zz = z; zz(~isfinite(zz)) = 0;
m = 0;
for s = [1 -1]
    above = (s*zz) > thr;                 % one sign at a time = sign-consistent runs
    if ~any(above), continue; end
    d = diff([false above false]);
    i1 = find(d == 1); i2 = find(d == -1) - 1;
    for k = 1:numel(i1)
        m = max(m, abs(sum(zz(i1(k):i2(k)))));
    end
end
end
