function batch_segment(subj_info_struct, base_dir)

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


for subj_num = 1:length(subj_info_struct) 
    this_subj_info = subj_info_struct(subj_num);
    fprintf('Running batch preprocessing on subject %d\n', ...
        this_subj_info.SubjectNum);
    %%% Load in the manually saved realign .mat file from subj 1's data
    load('lin_segment_template.mat');
    
    this_subj_base_dir = fullfile(imaging_data_subfolder, ...
        this_subj_info.ScanningDir, this_subj_info.ScanningDate, ...
        converted_dir );
    
    
    %%% Find all bold dirs for current subject (directories that contain
    %%% 'ep2d_bold' in their name and have at least 135 (num_TRs) files in
    %%% them
    this_subj_full_path = fullfile(base_dir, this_subj_base_dir);
    if ~ exist(this_subj_full_path, 'dir')
        fprintf('Converted directory not found: %s\n', this_subj_full_path);
        continue;
    end
    fprintf('Looking in %s for T1 directories \n', this_subj_full_path);
    subj_T1_dirs = findBoldDirs(this_subj_full_path, 1, 'T1_MPRAGE');
    
    T1_dir = fullfile(base_dir, this_subj_base_dir, subj_T1_dirs{1});
    %%% SPM select .nii files 
    [T1_image,dirs]=spm_select('ExtFPList', T1_dir, ...
                                       ['^' this_subj_info.ScanningDir '.*\.nii$']);
    matlabbatch{1}.spm.spatial.preproc.data{1} = T1_image;
    %%% Run the realign job
    
    spm_jobman('run',matlabbatch);
       
end;  % End of loop through subjects
