function multisvm_group_level(mode, type, vox)
% analyze group-level svm results

load(fullfile(BaseDir(), 'cache', 'classifier', sprintf('multisvm_%s_%s_stable%d.mat', ...
    type, mode, vox)));

 load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
    lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
    roi_table = struct2table(roi_masks_struct);
    names = roi_table.name;

if strcmp(mode, 'harvard-oxford')
    acc = zeros(96, length(sessions));
else
    acc = zeros(1, length(sessions));
end
for subj = 1:length(sessions)
    load(sessions(subj).classifier.mat);
    acc(:, subj) = mean(cv_acc);
end

% do t-test on each ROI, and finally do FDR correction
p_values = zeros(96, 1);
t_values = zeros(96, 1);
for roi = 1:96
    if strcmp(type, 'category')
        [h, p_values(roi), ci, stats] = ttest(acc(roi, :), 1/3, 'tail', 'Right');
    else
        [h, p_values(roi), ci, stats] = ttest(acc(roi, :), 1/15, 'tail', 'Right');
    end
    t_values(roi) = stats.tstat;
end

% do FDR correction
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values);
fprintf('left_brain: \n')
names(find(adj_p(1:48) < 0.05))
fprintf('right_brain: \n')
names(find(adj_p(49:96) < 0.05))
min_t = min(t_values(find(adj_p < 0.05)));
fprintf('minimum t-value is %.3f\n', min_t);

% write brain maps
% t-map
fname = fullfile(BaseDir(), 'figures', sprintf('Tmap_multisvm_%s_%s_stable%d.nii', ...
     type, mode, vox));
writeMask = fill_mask_voxels_with_roi_values(lr_rois, t_values);
write_brain(writeMask, lr_rois{1}, [], fname);
a=1;

