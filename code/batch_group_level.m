function [] = batch_group_level(subjects)
% batch_group_level(subjects)
%
% Run group-level analysis for <subjects> using currently active model.  Before
% running this, run batch_contrasts to make sure that the contrast files are
% generated.

if nargin < 1
    subjects = subj_info();
    warning('No subjects specified, using %d subjects from default list', length(subjects));
end

% load in SPM.mat files for each subject to extract the contrasts
subj_spms_cell = cell(size(subjects));
spm_filename_template = fullfile(...
    animalImagingBaseDir(), ...
    modelSubdir(), ...
    'new_subjects', ...
    'Lin%0.3d', ...
    '%s', ...
    'SPM_analysis', ...
    'SPM.mat');

for i = 1:length(subjects)
    subj_spms_cell{i} = load(sprintf(spm_filename_template, ...
        subjects(i).SubjectNum, subjects(i).PrePost));
    subj_spms_cell{i} = subj_spms_cell{i}.SPM;
end

% make struct array of SPMs
subj_spms = [subj_spms_cell{:}];
% make (subject x contrast) array of contrast structs.
subj_xcons = vertcat(subj_spms.xCon);
% make cell array of SPM working directories (where contrast images are stored)
subj_working_dirs = {subj_spms.swd};

contrast_names = {subj_xcons(1,:).name};
% check that every subject has all the contrasts in the right order
for i = 2:length(subjects)
    if ~all(strcmp(contrast_names, {subj_xcons(i,:).name}))
        error('Mismatch in contrast names for subject Lin%0.3d, %s. Try regenerating with batch_contrasts().', ...
            subjects(i).SubjectNum, subjects(i).PrePost);
    end
end

fprintf('Running group-level models on %d contrasts: \n', length(contrast_names));
fprintf('  %s\n', contrast_names{:});

% loop over contrasts
for contrast_idx = 1:size(subj_xcons,2)
    
    contrast_name = contrast_names{contrast_idx};
    
    % get list of absolute paths (in SPM file select format, so trailing ",1")
    % to contrast images for each subject
    contrast_image_filname = sprintf('con_%04d.img,1', contrast_idx);
    subject_filenames = fullfile(subj_working_dirs, contrast_image_filname);    
    
    % output directory: contrast name, sanitized
    illegal_chars = '[/\\"'';:()\s]';
    % replace any run of 1+ illegal characters with '_', and convert to
    % lowercase
    output_dir = lower(regexprep(contrast_name, [illegal_chars '+'], '_'));
    output_dir_full = fullfile(animalImagingBaseDir(), modelSubdir(), 'group_level', output_dir);
    
    % check for existence of output directory, create if not found
    if ~ exist(output_dir_full, 'dir')
        fprintf('Output directory not found\n  %s\n  creating...', output_dir_full);
        [success, message] = mkdir(output_dir_full);
        if ~ success
            error('Error creating output directory %s\n(%s)\n', output_dir_full, message);
        else
            fprintf('success!\n');
        end
    end
    
    config = struct();
    config.contrast_name = contrast_names{contrast_idx};
    config.contrast_images = subject_filenames;
    config.output_dir = output_dir_full;
    
    batch_obj = job_group_level(config);
    save(fullfile(output_dir_full, 'batch_group_level.mat'), 'batch_obj');
    
    % run the job
    try 
        fprintf('\n\n================================================================================\n');
        fprintf('Running group level analysis on contrast %s\n', contrast_name);
        spm_jobman('initcfg');
        spm_jobman('serial', batch_obj);
    catch exception
        fprintf('Encountered a problem:\n%s\n  Skipping this subject', ...
            getReport(exception));
    end
 

    
end