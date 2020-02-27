function batch_sliceTimeCorrect(subj_info_struct, base_dir)

%%% Find absolute path of base directory automagically

%%% Subdirectory of base_dir which contains imaging data. DON'T EDIT THIS
%%% if your data lives somewhere else; just create a symlink in base_dir
%%% which points to your data.  See README.md.
imaging_data_subfolder = 'data-imaging';

%%% Load default subject info if not provided
if nargin < 1
    [subj_info_struct, subj_info_cell] = subj_info();
end
if nargin < 2
    base_dir = BaseDir();
end
subj_dirs = {subj_info_struct.ScanningDir; 
             subj_info_struct.ScanningDate}';

% subj_dirs = { ...
%     'Lin002' '20131022'; ...  
%     'Lin002_post' '20131105'; ...
%     'Lin003' '20131112'; ...
%     'Lin004' '20131119'; ...   
%     'Lin004' '20131205'; ...
%     'Lin005' '20131212'
%    };

%%% Where do the NiFTi files live in each subject's directory?
converted_dir = 'Mcverter_Dicom_conversion/';

%%% How many TRs are there (used to filter out incomplete runs)
num_TRs = 170;


for subj_num = 1:length(subj_info_struct) 
    this_subj_info = subj_info_struct(subj_num);
    fprintf('Running batch preprocessing on subject %d\n', ...
        this_subj_info.SubjectNum);
    %%% Load in the manually saved normalize .mat file from subj 1's data
    load('lin_sliceTimeCorrect_template.mat');
    
    this_subj_base_dir = fullfile(imaging_data_subfolder, ...
        this_subj_info.ScanningDir, this_subj_info.ScanningDate, ...
        converted_dir );
    
    %%% Initialize list of files to empty cell array (will grow with each
    %%% bold dir processed)
    matlabbatch{1}.spm.temporal.st.scans = cell(1);
    
    %%% Find all bold dirs for current subject (directories that contain
    %%% 'ep2d_bold' in their name and have at least 135 (num_TRs) files in
    %%% them
    this_subj_full_path = fullfile(base_dir, this_subj_base_dir);
    if ~ exist(this_subj_full_path, 'dir')
        fprintf('Converted directory not found: %s\n', this_subj_full_path);
        continue;
    end
    fprintf('Looking in %s for BOLD directories \n', this_subj_full_path);
    subj_bold_dirs = findBoldDirs(this_subj_full_path, num_TRs);
    
    fprintf('Found %d BOLD directories: \n', length(subj_bold_dirs));

    %%% If there are no files to normalize, then should not run batch, so
    %%% keep track of whether any files are found
    any_file_to_normalize = false;
    
    %%% Iterate over directories, adding files to the job queue
    for this_bold_dir_idx = 1:length(subj_bold_dirs),
        fprintf('  BOLD Run %03d: %s\n', this_bold_dir_idx, subj_bold_dirs{this_bold_dir_idx});
        %%% Construct absolute path for currently processed directory
        this_bold_dir_string = fullfile(base_dir, this_subj_base_dir, subj_bold_dirs{this_bold_dir_idx});
        %%% SPM select .nii files 
        [these_bold_files,dirs]=spm_select('ExtFPList',this_bold_dir_string, ...
                                           ['^r' this_subj_info.ScanningDir '.*\.nii$']);
        %%% Store in cell array to add to jobs struct
        bold_files_cell = mat2cell(these_bold_files,ones(1,size(these_bold_files, 1)),size(these_bold_files,2));
        
        %%% Check to see if normalized files are already present.
        [paths, fns] = cellfun(@fileparts, bold_files_cell, 'UniformOutput', false);
        normalized_bold_files_exist = cellfun(@(path, fn) ...
            logical(exist(fullfile(path, ['a' fn '.nii']))), paths, fns);
        if any(normalized_bold_files_exist)
            fprintf('    %d / %d .nii files have normalized versions already!\n', ...
                    sum(normalized_bold_files_exist), length(normalized_bold_files_exist));
        end

        %%% Put this subject's BOLD files (not already aligned) into the jobs struct
        if ~ all(normalized_bold_files_exist)
            matlabbatch{1}.spm.temporal.st.scans{1} = ...
                [matlabbatch{1}.spm.temporal.st.scans{1}; ...
                 bold_files_cell(~normalized_bold_files_exist)];
            any_file_to_normalize = true;
        else
            fprintf('    No files to re-align.  Delete existing files to re-align again.\n');
        end
    end; % end bold dir loop
    
    
    %%% Run the normalize job
    if any_file_to_normalize
        spm_jobman('run',matlabbatch);
    end
       
end;  % End of loop through subjects
