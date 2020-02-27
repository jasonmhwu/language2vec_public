% setup script

% enter path you store spm12 and 
spm12_path = '/home/mwu34/Documents/MATLAB';
addpath(fullfile(spm12_path));


% load all relevant folders
root_path = fileparts(fileparts(mfilename('fullpath'))); % base directory
addpath(genpath(fullfile(root_path, 'code')));
addpath(genpath(fullfile(root_path, 'public_repositories')));
