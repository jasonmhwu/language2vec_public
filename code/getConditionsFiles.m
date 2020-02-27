function [conditionsFiles, subj_info] = getConditionsFiles(subj_info, generate_if_missing)
% [conditionsFiles, subj_info] = getConditionsFiles(subj_info, generate_if_missing)
%
% return cell array of conditions files for a given subject (or structure
% array of subjects).  Optionally returns an updated subject info struct
% with field `ConditionsFiles` added.

if nargin < 2
    generate_if_missing = 0;
end

conditions_file_dir_template = fullfile(BaseDir(), ...
                                        'subjects', ...
                                        'Subject%03d', ...
                                        'conditions');
matfile_name_template = 'conditions_run%03d.mat';



for subject_idx = 1:length(subj_info)
       
    this_subj_info = subj_info(subject_idx);
    
    % get runs 
    runs = parseOneSession(this_subj_info.SubjectNum);
    nruns = length(runs);
    
    % initialize cell array of conditions file names
    subjConditionsFiles = cell(1, nruns);
    
    conditions_file_dir = sprintf(conditions_file_dir_template, ...
                                  this_subj_info.SubjectNum);
    if ~exist(conditions_file_dir)
        mkdir(conditions_file_dir);
    end
    % iterate over runs
    for run_num = 1:nruns
        run_matfile_name = sprintf(matfile_name_template, run_num);
        conditions_matfile_absolute = fullfile(conditions_file_dir, run_matfile_name);

        % generate conditions files if missing (will use model-specific
        % makeConditionsFiles() function)
        if generate_if_missing && ~ exist(conditions_matfile_absolute)
            fprintf('WARNING: conditions file does not exist.  Will try to generate.\n');
            onsets = num2cell(cell2mat(runs{run_num}.data(strcmp(runs{run_num}.varnames, 'onsetTime'))));
            durations = num2cell(cell2mat(runs{run_num}.data(strcmp(runs{run_num}.varnames, 'duration'))));
            for i = 1:length(onsets)
                names{i} = sprintf('run-%d-trial-%d', run_num, i);
            end
            names = reshape(names, length(names), 1);
            save(conditions_matfile_absolute, 'onsets', 'durations', 'names');
        end
        
        subjConditionsFiles{run_num} = conditions_matfile_absolute;
    end
    
    conditionsFiles{subject_idx} = subjConditionsFiles;
    subj_info(subject_idx).ConditionsFiles = subjConditionsFiles;
end