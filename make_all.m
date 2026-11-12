function make_all(varargin)
% make_all  Rebuild Figures 7 and 8 (and the two supplementary videos) from data/.
%
%   make_all                    figures only (a few minutes)
%   make_all('Videos', true)    figures and videos
%
% Output goes to output/. Run verify afterwards to check the numbers.
p = inputParser;
p.addParameter('Videos', false);
p.parse(varargin{:});
setup_paths();
st = jneStyle();
D = jneData();
for f = {@fig7c_polarity_map, @fig7d_grid_waveforms, @fig7e_p2p_progression, ...
         @fig7f_psd_progression, @fig7g_channel_yield}
    fprintf('[make_all] %s\n', func2str(f{1}));
    close(f{1}(D, st));
end
jne_contactInset(D, st); close all;
fprintf('[make_all] fig8_travelingWaves\n');
close(fig8_travelingWaves());
if p.Results.Videos, jne_suppVideos(); end
fprintf('[make_all] done -> %s\n', jneOutDir());
end
