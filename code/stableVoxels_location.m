function stableVoxels_location(sessions, Nvox)
% estimate stable Voxels for each wordType
% should do standardization first

if nargin < 1
    sessions = subj_info();
end

if nargin < 2
    Nvox = 800;
end

output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'stableVoxels');

% load harvard-oxford atlas
load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
lr_rois = {roi_masks_struct.left, roi_masks_struct.right};

stableVoxCount = zeros(length(lr_rois), length(sessions));
actualVoxCount = zeros(length(lr_rois), length(sessions));

for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    stable_fname = sprintf(fullfile(output_dir_template, 'corrXRuns.nii'), ...
        session.SubjectNum);
    stableVoxels = spm_read_vols(spm_vol(stable_fname));
    [all_betas_mat, runInfo, mask_mat] = ...
            retrieve_actualTrials_one_sub(session);
    
    
     stabilityInRoi =sort(stableVoxels(logical(mask_mat)), 'descend');
     stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
     
     for roi_idx = 1:length(lr_rois)
     stableVoxCount(roi_idx, sess_idx) = length(find(lr_rois{roi_idx} & ...
         (stableVoxels >= stabilityInRoi(min(Nvox, length(stabilityInRoi))))));
     actualVoxCount(roi_idx, sess_idx) = length(find( ...
         lr_rois{roi_idx}& logical(mask_mat)));
     end
end

plot(mean(stableVoxCount, 2));
[~, mostStableROI] = sort(mean(stableVoxCount, 2), 'descend');

% divided by number of voxels in each ROI?
[~, mostStableROI_normalized] = sort( ...
    mean(stableVoxCount ./ actualVoxCount, 2), 'descend');
a=1;
end