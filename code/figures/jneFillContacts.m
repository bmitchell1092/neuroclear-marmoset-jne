function vals = jneFillContacts(vals, chMap, drop)
% jneFillContacts(vals, chMap, drop)
% Replace dropped contacts' values with a natural-neighbour interpolation from
% the surviving contacts (at their array grid positions), keeping ALL contact
% positions finite so the array footprint is unchanged. 'nearest' extrapolation
% covers dropped contacts that sit outside the good-contact hull (array edges).
% Shared by fig7c_polarity_map and fig7e_p2p_progression.
%
%   vals  : [nChan x 1] per-contact metric (erp-row order)
%   chMap : channelMap (0-based channel index per grid cell; NaN = no contact)
%   drop  : [nChan x 1] logical, true = fill this contact from neighbours
[nr, nc] = size(chMap);
gr=[]; gc=[]; gv=[]; br=[]; bc=[]; bidx=[];
for rr = 1:nr
    for cc = 1:nc
        ch = chMap(rr,cc); if isnan(ch), continue; end
        e = ch + 1; if e < 1 || e > numel(vals), continue; end
        if drop(e) || ~isfinite(vals(e))
            br(end+1)=rr; bc(end+1)=cc; bidx(end+1)=e; %#ok<AGROW>
        else
            gr(end+1)=rr; gc(end+1)=cc; gv(end+1)=vals(e); %#ok<AGROW>
        end
    end
end
if numel(gv) >= 3 && ~isempty(bidx)
    F = scatteredInterpolant(gc(:), gr(:), gv(:), 'natural', 'nearest');
    vals(bidx) = F(bc(:), br(:));
end
end
