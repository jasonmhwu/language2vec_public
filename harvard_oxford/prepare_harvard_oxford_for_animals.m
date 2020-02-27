% depends on raizadalab-utilities: 
% https://bitbucket.org/raizadalabteam/raizadalab-utilities

% NOTE: the outputs are included in the repo for maximum portability, but this
% script is here just to document how they were generated.

%% create realigned atlas image
ho_fn = fullfile(BaseDir(), ...
                 'harvard_oxford', ...
                 'HarvardOxford-cort-maxprob-thr0-2mm.nii');
target_beta_fn = fullfile(BaseDir(),...
                          'subjects/Subject001/SPM_analysis/beta_0001.nii');
output_fn = fullfile(BaseDir(), ...
                     'harvard_oxford', ...
                     'harvard_oxford_animal_space.nii');

put_into_same_voxel_space_spm5(ho_fn, target_beta_fn, output_fn);

%% clean up and geneate ROI images/struct

harv_oxf_dir = fullfile(BaseDir(), 'harvard_oxford');
% strip dir off output filename
[~, harv_oxf_realigned_fn] = fileparts(output_fn);

make_harvard_oxford_roi_masks(harv_oxf_dir, [harv_oxf_realigned_fn '.nii']);

%% create realigned atlas image - subcortical
ho_fn = fullfile(BaseDir(), ...
                 'harvard_oxford', ...
                 'HarvardOxford-sub-maxprob-thr0-2mm.nii');
target_beta_fn = fullfile(BaseDir(),...
                          'subjects/Subject001/SPM_analysis/beta_0001.nii');
output_fn = fullfile(BaseDir(), ...
                     'harvard_oxford', ...
                     'harvard_oxford_sub_animal_space.nii');

put_into_same_voxel_space_spm5(ho_fn, target_beta_fn, output_fn);
% write to a mat file
rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
                                     'harvard_oxford_sub_animal_space.nii', '', 0);
                                 