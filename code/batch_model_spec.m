function batch_model_spec(subj_info_struct, conditions_mat_files_subdir)

%%% Find absolute path of base directory automagically
base_dir = animalImagingBaseDir();
%%% Subdirectory of base_dir which contains imaging data. DON'T EDIT THIS
%%% if your data lives somewhere else; just create a symlink in base_dir
%%% which points to your data.  See README.md.
imaging_data_subfolder = 'data-imaging';

%%% Load default subject information if not provided.
if nargin < 1
    subj_info_struct = subj_info();
end

subj_dirs = {subj_info_struct.ScanningDir; 
             subj_info_struct.ScanningDate}';

%%% Where do the NiFTi files live in each subject's directory?
converted_subdir = 'new_Mcverter_Dicom_conversion/';
%%% What is the regex pattern to select data files? 
image_file_regex = '^war.*\.nii$';

%%% Where does the SPM analysis happen? This is where teh SPM.mat file will
%%% be written, and is a subdirectory of converted_subdir
analysis_subdir = 'SPM_analysis_batch/';

%%% where are the conditions .mat files stored? this is a sub directory of
%%% the project base directory.
if nargin < 2
    conditions_mat_files_subdir = 'conditions-files/';
end

%%% How many TRs are there (used to filter out incomplete runs)
num_TRs = 135;

for subj_num = 1:length(subj_info_struct), 
    fprintf('Running batch preprocessing on subject %d\n', subj_num);
    %%% Load in the manually saved fmri_spec .mat file from subj 1's data
    load('lin_fMRI_spec_template.mat');
    
    this_subj_info = subj_info_struct(subj_num);
    
    %%% cd to subject directory
    this_subj_base_dir = fullfile(imaging_data_subfolder, ...
        this_subj_info.ScanningDir, this_subj_info.ScanningDate, converted_subdir );

    %%% Find all bold dirs for current subject (directories that contain
    %%% 'ep2d_bold' in their name and have at least 135 (num_TRs) files in
    %%% them
    this_subj_full_path = fullfile(base_dir, this_subj_base_dir);
    fprintf('Looking in %s for BOLD directories \n', this_subj_full_path);
    subj_bold_dirs = findBoldDirs(this_subj_full_path, num_TRs);
    
    if length(subj_bold_dirs) > 0

        fprintf('Found %d BOLD directories: \n', length(subj_bold_dirs));

        %%% The only thing we need to change is that we are going to
        %%% load in a different subject's data
        bold_files_from_manual_job = matlabbatch{1}.spm.stats.fmri_spec.sess;  
        all_bold_files = cell(1,num_TRs);
        for this_bold_dir_idx = 1:length(subj_bold_dirs),
            fprintf('  BOLD Run %03d: %s\n', this_bold_dir_idx, subj_bold_dirs{this_bold_dir_idx});
            %%% Construct absolute path for currently processed directory
            this_bold_dir_string = fullfile(base_dir, this_subj_base_dir, subj_bold_dirs{this_bold_dir_idx});
            %%% SPM select realigned .nii files (begin with 'r'). 
            %%% In fact, we will actually want the normalised ones, starting
            %%% with nr, and the smoothed normalised snr
            [these_bold_files,dirs]=spm_select('ExtFPList',this_bold_dir_string, ...
                                               image_file_regex);
            bold_files_cell = mat2cell(these_bold_files,ones(1,num_TRs),size(these_bold_files,2));

            %%% Put this subject's BOLD files into the jobs struct
            matlabbatch{1}.spm.stats.fmri_spec.sess(this_bold_dir_idx).scans = bold_files_cell;

            %%% Load in conditions .mat file

            %%% Get the conditions .mat files, 
            %%% Auto-generate conditions files if they don't exist.
            onsets_matfile_path = fullfile(base_dir, conditions_mat_files_subdir, ...
                sprintf('Lin%03d_%s', this_subj_info.SubjectNum, this_subj_info.PrePost));
            matfile_name = sprintf('conditions_run%d.mat', this_bold_dir_idx);
            conditions_matfile_absolute = fullfile(onsets_matfile_path, matfile_name);
            if ~ exist(conditions_matfile_absolute)
                fprintf('WARNING: conditions file does not exist.  Will try to generate.\n');
                makeConditionsFiles(this_subj_info.SubjectNum, this_subj_info.PrePost, output_dir_template);
            end

            %%% Record file name in batch structure
            matlabbatch{1}.spm.stats.fmri_spec.sess(this_bold_dir_idx).multi{1} = ...
                conditions_matfile_absolute;        

        end;  % End of loop through runs

        %%% Set the directory where the the SPM.mat file will be written
        %%% (creating it first if necessary)
        analysis_dir_absolute = fullfile(this_subj_full_path, analysis_subdir);
        if ~ exist(analysis_dir_absolute, 'dir')
            mkdir(analysis_dir_absolute);
        end
        matlabbatch{1}.spm.stats.fmri_spec.dir{1} = ...
            fullfile(this_subj_full_path, analysis_subdir);

        %%% Run the fMRI_spec job
        spm_jobman('run',matlabbatch);
    else
        %%% No BOLD directories found, abort
        fprintf('No bold directories found, skipping.\n');
    end
        
end;  % End of loop through subjects
