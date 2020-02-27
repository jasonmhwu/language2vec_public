function dirname = BaseDir()

% Returns the directory containing this m file. Short-cut to get the base
% directory for the analysis code.
[dirname, filename] = fileparts(fileparts(mfilename('fullpath')));
