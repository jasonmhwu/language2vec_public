function boldDirs = findBoldDirs(parentDir, minFileCount, boldPattern)
% function boldDirs = findBoldDirs(parentDir, boldPattern, minFileCount)
%
% find all the subdirectories of the specified parent directory that are
% BOLD runs.
% 
% Optional arguments boldPattern='ep2d_bold' by default, and minFileCount=0
% is the number of files in the directory for it to be returned

% default pattern
if nargin < 3
    boldPattern = 'ep2d_bold';
end

% make is in this directory
contents = dir(parentDir);

% return names of all directories 
dirNames = {contents([contents(:).isdir]).name};

% does the folder name match the boldPattern?
boldDirInds = cellfun(@(x) (~isempty(x)), ...
                      regexp(dirNames, boldPattern));
                  
% produces a list of folder names matching the boldPattern
boldDirs = dirNames(boldDirInds);

% produces a list of boldDirs that contains more than minFileCount files (+
% 2 files for '.' and '..')
if nargin > 1 && minFileCount > 0
    enoughContentsInds = ...
        cellfun(@(x) (length(dir(fullfile(parentDir, x))) > minFileCount+2), boldDirs);
    boldDirs = boldDirs(enoughContentsInds);
end