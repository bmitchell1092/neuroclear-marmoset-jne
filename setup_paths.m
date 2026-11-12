function setup_paths()
% setup_paths  Put this repo's code and the bundled toolboxes on the MATLAB path.
r = fileparts(mfilename('fullpath'));
ext = fullfile(r, 'external');
if ~isfile(fullfile(ext, 'wave-matlab', 'analysis', 'analytic_signal.m'))
    error('jne:submodules', ['external/wave-matlab is empty. Clone with --recursive, or run\n' ...
        '  git submodule update --init']);
end
addpath(fullfile(r,'code','data'), fullfile(r,'code','figures'), fullfile(r,'code','wave'), ...
    fullfile(ext,'wave-matlab','analysis'), fullfile(ext,'generalized-phase','analysis'), ...
    fullfile(ext,'smoothn'), fullfile(ext,'circstat'));
end
