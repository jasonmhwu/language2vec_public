function writeMask = fill_mask_voxels_with_roi_values(rois, values, mask)
% Map values onto voxels according to ROIs
% 
% [voxel_vec] = fill_mask_voxels_with_roi_values(rois, values, mask)
% 
% Inputs: 
%   rois: cell array of ROIs in mask indices
%   values: vector of value for each ROI
%   mask: binary (0/1) mask image, or number of voxels contained in mask. (only
%     used to initialize the output vector to the right length).
% Outputs: 
%   voxel_vec: A vector of values for each voxel in the mask.
% 
% Example:
% 
% >> roi_acc_map = fill_mask_voxels_with_roi_values(rois, roi_acc, mask);
% >> write_brain(roi_acc_map, mask, header, 'my_roi_acc_map.nii');
if nargin < 3
    mask = spm_read_vols(spm_vol(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space.nii')));
end

if length(values) ~= length(rois)
    error('Size mismatch between values and number of rois');
end

writeMask = nan(size(mask));

for roi_idx = 1:length(rois)
    writeMask(rois{roi_idx}) = values(roi_idx);
end

writeMask(mask == 0) = nan;