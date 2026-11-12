function featCh = jneResolveContact(sel, D, chMap, n)
% jneResolveContact(sel, D, chMap, n) -> erp-row index (1-based) of a contact.
% Accepts:
%   []            auto: contact with the strongest per-session-minimum p2p
%   scalar        an erp-row index (1..nChan), returned as-is
%   [row col]     a grid location in channelMap (0-based channel -> erp row+1)
% Shared by fig7d_grid_waveforms and jne_contactInset so "which contact" is
% specified the same way everywhere.
if nargin < 4 || isempty(n), n = numel(D.flash); end
if isempty(sel)
    M = nan(D.flash(1).R.nChan, n);
    for i = 1:n, M(:,i) = D.flash(i).R.ampP2P(:); end
    [~, featCh] = max(min(M, [], 2));
elseif numel(sel) == 2
    ch0 = chMap(sel(1), sel(2));
    if isnan(ch0), error('jneResolveContact:badContact', ...
            'No contact at grid (row %d, col %d).', sel(1), sel(2)); end
    featCh = ch0 + 1;
else
    featCh = sel;
end
end
