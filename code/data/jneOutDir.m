function d = jneOutDir()
% jneOutDir  output/ at the top of the repo (created if missing)
d = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'output');
if ~isfolder(d), mkdir(d); end
end
