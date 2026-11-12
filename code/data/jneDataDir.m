function d = jneDataDir()
% jneDataDir  data/ at the top of the repo (download it from Zenodo, see README)
d = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data');
if ~isfolder(fullfile(d, 'fig7')) || ~isfolder(fullfile(d, 'fig8'))
    error('jne:noData', 'No data in %s. Download it from Zenodo and unzip it there (see README).', d);
end
end
