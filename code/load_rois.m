function [lr_rois, names] = load_rois(mode)
    names = {};

    if strcmp(mode, 'harvard-oxford')
        %rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
        %                                 'harvard_oxford_animal_space.nii');
        load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
        lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
        roi_table = struct2table(roi_masks_struct);
        names = roi_table.name;
    elseif strcmp(mode, 'fedorenko-rois')
        load(fullfile(BaseDir(), 'fedorenko_rois', 'fedorenko_rois_realigned_all_rois.mat'));
        lr_rois = {roi_masks_struct.left};
        roi_table = struct2table(roi_masks_struct);
        names = roi_table.name(rois);
    elseif strcmp(mode, 'GroupMask_stable200_union8')
        load(fullfile(BaseDir(), 'cache', 'GroupMask_word2vec_fedorenko-rois_stable200_intersect8.mat'));
        lr_rois = GroupMask_union;
    elseif strcmp(mode, 'union-language-network')
        load(fullfile(BaseDir(), 'fedorenko_rois', 'fedorenko_rois_realigned_all_rois.mat'));
        union_roi = zeros(size(roi_masks_struct(1).left));
        for roi = 1:length(roi_masks_struct)
            union_roi = union_roi + double(roi_masks_struct(roi).left);
        end
        union_roi = logical(union_roi);
        lr_rois{1} = union_roi;
        lr_rois{2} = ~union_roi;
        names = {'inNetwork', 'outNetwork'};
    elseif contains(mode, 'general-semantic-network')
        % OPTION: load from HO
        %{ 
        load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
        GSN_idx = {[19, 20, 21], ...
                           [11, 12, 13, 15], ...
                           [34, 35, 37, 38, 39], ...
                           [3, 4], ...
                           [5, 6], ...
                           [1, 33], ...
                           [30, 31], ...
                           [8]};
        %}
        
        % OPTION: load from AAL
        load(fullfile(BaseDir(), 'aal', 'aal_all_rois.mat'));
        GSN_idx = {[34, 35], ...
                            [45, 47], ...
                            [30, 22], ...
                            [2, 10], ...
                            [4, 5, 6], ...
                            [11, 12, 13], ...
                            [20, 36], ...
                            [44, 46]};
        for roi = 1:length(GSN_idx)
            union_roi = zeros(size(roi_masks_struct(1).left));
            for r = GSN_idx{roi}
                union_roi = union_roi + double(roi_masks_struct(r).left);
            end
            lr_rois{roi} = logical(union_roi);
        end
        % option to write out map
        %{
        fullROI = zeros(size(roi_masks_struct(1).left));
        for roi = 1:length(GSN_idx)
            fullROI = fullROI + roi * double(lr_rois{roi});
        end
        write_brain(fullROI, fullROI, [], ...
            fullfile(BaseDir(), 'aal', 'binder_GSN_aal.nii'));
        %}
        names = {'IPL', 'LTL', 'VTL', 'dmPFC', 'IFG', 'vmPFC', 'PCG', 'ATL'};
    elseif strcmp(mode, 'harvard-oxford-cortical-subcortical');
        load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
        union_roi = zeros(size(roi_masks_struct(1).left));
        for roi = 1:length(roi_masks_struct)
            union_roi = union_roi + double(roi_masks_struct(roi).left) + double(roi_masks_struct(roi).right);
        end
        lr_rois{1} = logical(union_roi);
        load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_sub_animal_space_all_rois.mat'));
        union_roi = zeros(size(roi_masks_struct(1).left));
        for roi = [1, 3, 4, 5, 6, 7, 9, 10]
        %for roi = 1
            union_roi = union_roi + double(roi_masks_struct(roi).left) + double(roi_masks_struct(roi).right);
        end
        lr_rois{2} = logical(union_roi);
    elseif contains(mode, 'MNI-four-lobes')
        load(fullfile(BaseDir(), 'MNI', 'MNI-space_all_rois.mat'));
        rois = [3, 5, 6, 8];
        lr_rois = {roi_masks_struct(rois).left, roi_masks_struct(rois).right};
        roi_table = struct2table(roi_masks_struct);
        names = roi_table.name(rois);
    else % whole-brain
        lr_rois = cell(1);
    end
end