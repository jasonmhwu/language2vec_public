% depends on raizadalab-utilities: 
% https://bitbucket.org/raizadalabteam/raizadalab-utilities

% NOTE: the outputs are included in the repo for maximum portability, but this
% script is here just to document how they were generated.

%% create realigned atlas image
atlas_dir = '/raizadaUsers/mwu34/Documents/MATLAB/spm12/toolbox/aal/atlas';
ho_fn = fullfile(atlas_dir, 'AAL2.nii');
target_beta_fn = fullfile(BaseDir(),...
                          'subjects/Subject001/SPM_analysis/beta_0001.nii');
output_dir = fullfile(BaseDir(), 'aal');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_fn = fullfile(BaseDir(), ...
                     'aal', ...
                     'aal.nii');

put_into_same_voxel_space_spm5(ho_fn, target_beta_fn, output_fn);

%% clean up and geneate ROI images/struct

% strip dir off output filename
[~, aal_realigned_fn] = fileparts(output_fn);

make_aal_roi_masks(output_dir, [aal_realigned_fn '.nii']);
