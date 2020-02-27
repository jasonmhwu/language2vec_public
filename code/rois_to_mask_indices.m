function [roi_mask_indices] = rois_to_mask_indices(rois, mask, replaceNaN)
% function [roi_mask_indices] = rois_to_mask_indices(rois, mask)
% 
% Convert ROIs in mask format (volumes with non-zero values indicating voxels in
% the ROI) to ROIs in mask-space index format (vector of indices in a different
% mask).  This allows ROIs to be used to select subsets of columns in a trial x
% voxel matrix of beta values such as that produced by assemble_betas_one_sub.

if nargin < 3
    replaceNaN = false;
end

if ~islogical(mask)
    mask = logical(mask);
end

roi_mask_indices = cell(size(rois));

for roi_num = 1:length(rois)
    this_roi = rois{roi_num};
    if ~replaceNaN
        roi_mask_indices{roi_num} = find(this_roi(mask));
    else
        InMask = find(this_roi(mask));
        NotInMask = find(this_roi(~mask));
        tmp_indices = sort([InMask; NotInMask]);
        for i = 1:length(tmp_indices)
            if ismember(tmp_indices(i), NotInMask)
                tmp_indices(i) = NaN;
            end
        end
        roi_mask_indices{roi_num} = tmp_indices;
        
    end
end
