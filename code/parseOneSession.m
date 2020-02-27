function dataRuns = parseOneSession(subjectnumber) 
% function dataRuns = parseOneSession(subjectnumber, prepost, blocks) 
%   subjectnumber: duh.
%   prepost: which session to parse ('pre' or 'post')
%   blocks: cell array of blocks to parse (generalization, pitch, and passive
%           by default
%
% For each block in blocks, loads 
%   LinLearning_<prepost>scan_<block>_data_<subjectnumber>.txt
% parses it, and separates out the runs.  Then it sorts all the blocks by
% actual start time, making a cell array of the individual runs with trial
% times relative to start of run (in seconds and---possibly
% fractional---TRs)


% get first trial time to sort runs
%firstTs = zeros(size(datas));

base_dir = BaseDir();
behavioral_data_subdir = 'behavioral/';

% fprintf('Parsing data files: \n');
foldername = ...
        fullfile(base_dir, behavioral_data_subdir, ...
            sprintf('subj_%02d', subjectnumber));

% find txt files, and sort according to timestamp
fileStruct = dir([foldername, '/*.txt']);
timeStamp = string();
for i = 1:length(fileStruct)
    C = textscan(fileStruct(1).name,'%*s %*s %*s %*s %s %s %s %s %s %s %*s','Delimiter','_');
    for j = 1:length(C)
        timeStamp(i) = [C{2}{1}, C{3}{1}, C{4}{1}, C{5}{1}, C{6}{1}];
    end
end
[~, idx] = sort(timeStamp);
fileStruct = fileStruct(idx);

% parse each file
for i=1:length(fileStruct)
    % create the file name
    
    % parse the data file, and put in the cell array of parsed data files
    dataRuns{i} = parseScanBehavData(fullfile(fileStruct(i).folder, fileStruct(i).name));
    dataRuns{i}.runs = i;
end
dataRuns = dataRuns';