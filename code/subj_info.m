function [subj_info_struct, subj_info_cell] = subj_info(loadCatchTrial)
% function [subj_info_struct, subj_info_cell] = subj_info()
%
% Load subject information.  The struct has fields: 
%   PrePost = 'pre' or 'post'
%   ScanningDate = '20140130' (date of scan, subdir of ScanningDir)
%   ScanningDir = 'Lin004' (directory where scan data is actually found)
%   SubjectNum = 4 (the real subject number)
% 
% Reads in the file subj_info.txt which is generated from
% extract_info_from_behav_data.sh which matches the date and subject number
% from the behavioral data files with the date and subject numbers from the
% imaging data subdirectories. This shell script has to be run manually
% before running this function. At the scanner some subject numbers were
% entered incorrectly so the date and the subject provides redundant
% information to help correct these kind of errors.
if nargin < 1
    loadCatchTrial = 0;
end
info = textscan(fopen(fullfile(BaseCodeDir(), 'subj_info.txt'), 'r'), ...
    '%d %s %s %d %d %d %d\n');
subj_info_cell = [num2cell(info{1}), info{2}, info{3}];
subj_info_colnames = {'SubjectNum', 'ScanningDir', 'ScanningDate'};

subj_info_struct = cell2struct(subj_info_cell, subj_info_colnames, 2);

if loadCatchTrial
    subj_info_cell = [num2cell(info{1}), info{2}, info{3}, ...
        num2cell(info{4}), num2cell(info{5}), num2cell(info{6}), num2cell(info{7})];
subj_info_colnames = {'SubjectNum', 'ScanningDir', 'ScanningDate', ...
    'Correct', 'Incorrect', 'Miss', 'Error'};

subj_info_struct = cell2struct(subj_info_cell, subj_info_colnames, 2);

end