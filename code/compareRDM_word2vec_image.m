function compareRDM_word2vec_image(sessions, mode, wordModel, NstableVox, Nperm)
% can simple addition be a good model for specific areas?
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

if nargin < 1
    sessions = subj_info();
end
if nargin < 2
    mode = 'harvard-oxford';
end
if nargin < 3
    wordModel = 'word2vec';
end
if nargin < 4
    NstableVox = 100;
end
if nargin < 5
    Nperm = 1000;
end

if strcmp(mode, 'harvard-oxford')
    %rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
    %                                 'harvard_oxford_animal_space.nii');
    load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
    lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
end

output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'compareRDM');
stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');
    
for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    [all_betas_mat, all_betas_info, mask_mat] = ...
            assemble_betas_one_sub(session.SubjectNum, 0, 0, 'SPM_analysis');
        
    output_dir = sprintf(output_dir_template, session.SubjectNum);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % standardize images
    runInfo = parseAndCombineOneSession(session.SubjectNum);
    runIdx = getColumn(runInfo, 'runNum');
    runIdx_uniq = unique(runIdx);
    for i = 1:length(runIdx_uniq)
        tmp = all_betas_mat(find(runIdx == runIdx_uniq(i)), :);
        tmp = (tmp - mean(tmp)) ./ std(tmp);
        all_betas_mat(find(runIdx == runIdx_uniq(i)), :) = tmp;
    end
    
    
    wordType = getColumn(runInfo, 'wordType');
    Stimulus = getColumn(runInfo, 'Stimulus');
    stimCategory = getColumn(runInfo, 'stimCategory');
    
    % find the most stable voxels
    stable_fname = sprintf(stableVoxel_fname_template, session.SubjectNum);
    stableVoxels = spm_read_vols(spm_vol(stable_fname));
    if strcmp(mode, 'harvard-oxford')
        stable_lr_rois = cell(1, length(lr_rois));
        for i = 1:length(lr_rois)
            stabilityInRoi =sort(stableVoxels(lr_rois{i}), 'descend');
            stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
            stable_lr_rois{i} = lr_rois{i} & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
        end
    else % store whole brain as one ROI
        stabilityInRoi =sort(stableVoxels(logical(mask_mat)), 'descend');
        stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
        stable_lr_rois{1} = mask_mat & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
    end
    % use harvard-oxford atlas
    roi_mask_indices = rois_to_mask_indices(stable_lr_rois, mask_mat);
    corrRDM = zeros(2, length(stable_lr_rois));
    p_value = zeros(2, length(stable_lr_rois));
    for roi_idx = 1:length(stable_lr_rois)
        beta_masked = all_betas_mat(:, roi_mask_indices{roi_idx});
    
        % average across runs to create beta matrix
        Stimulus_uniq = unique(Stimulus);
        stimMatrix = struct();
        for i = 1:length(Stimulus_uniq)
            trials = strcmp(Stimulus, Stimulus_uniq{i});
            stimMatrix(i).beta = mean(beta_masked(trials, :));
            stimMatrix(i).wordType = wordType(find(trials, 1));
            stimMatrix(i).stimCategory = stimCategory(find(trials, 1));
            stimMatrix(i).Stimulus = Stimulus_uniq{i};
        end
        stimMatrix = struct2table(stimMatrix);

        % load wordModel
        [RDM, stimMatrix] = load_wordModel_RDM(stimMatrix, wordModel);
        %visualizeMDS(stimMatrix, 'wordModel_vec', fullfile(output_dir, 'word2vec_mds.jpg'));
        
        % calculate real fMRI RDM
        roi_RDM = corr(stimMatrix.beta');
        if strcmp(mode, 'harvard-oxford')
            if roi_idx <= 48
                visualizeMDS(stimMatrix, 'beta', sprintf(fullfile(output_dir, 'roi_left_%s_mds.jpg'), ...
                    roi_masks_struct(roi_idx).name));
            else
                visualizeMDS(stimMatrix, 'beta', sprintf(fullfile(output_dir, 'roi_right_%s_mds.jpg'), ...
                    roi_masks_struct(roi_idx-48).name));
            end
        else
            visualizeMDS(stimMatrix, 'beta', fullfile(output_dir, 'whole_brain_mds.jpg'));
        end
        % calculate Pearson & Spearman correlation
        corrRDM(1, roi_idx) = corr(roi_RDM(triu(true(size(roi_RDM)),1)), ...
                                   RDM(triu(true(size(RDM)),1)), ...
                                   'type', 'Pearson');
                               
        corrRDM(2, roi_idx) = corr(roi_RDM(triu(true(size(roi_RDM)),1)), ...
                                   RDM(triu(true(size(RDM)),1)), ...
                                   'type', 'Spearman');
        
        % permutation test
        permCorr = zeros(2, Nperm);
        for perm = 1:Nperm
            permIdx = randperm(size(RDM, 1));
            perm_RDM = RDM(permIdx, permIdx);
            permCorr(1, perm) = corr(roi_RDM(triu(true(size(roi_RDM)),1)), ...
                                   perm_RDM(triu(true(size(perm_RDM)),1)), ...
                                   'type', 'Pearson');
            permCorr(2, perm) = corr(roi_RDM(triu(true(size(roi_RDM)),1)), ...
                                   perm_RDM(triu(true(size(perm_RDM)),1)), ...
                                   'type', 'Spearman');
        end
        p_value(1, roi_idx) = sum(permCorr(1, :) > corrRDM(1, roi_idx)) / Nperm;
        p_value(2, roi_idx) = sum(permCorr(2, :) > corrRDM(2, roi_idx)) / Nperm;
        % report results
        if strcmp(mode, 'harvard-oxford')
            if roi_idx <= 48
                fprintf('left %s has %d voxels.\n', roi_masks_struct(roi_idx).name, size(beta_masked, 2));
            else
                fprintf('right %s has %d voxels.\n', roi_masks_struct(roi_idx-48).name, size(beta_masked, 2));
            end
            fprintf('%s model have an Pearson Correlation of %f,  permutation test p-value is %f\n', ...
                wordModel, corrRDM(1, roi_idx), p_value(1, roi_idx));
            fprintf('%s model have an Spearman Correlation of %f,  permutation test p-value is %f\n', ...
                wordModel, corrRDM(2, roi_idx), p_value(2, roi_idx));
        else
            fprintf('whole brain has %d voxels.\n', size(beta_masked, 2));
            fprintf('%s model have an Pearson Correlation of %f,  permutation test p-value is %f\n', ...
                wordModel, corrRDM(1, roi_idx), p_value(1, roi_idx));
            fprintf('%s model have an Spearman Correlation of %f,  permutation test p-value is %f\n', ...
                wordModel, corrRDM(2, roi_idx), p_value(2, roi_idx));
        end
    end
    
    
    
    % fdr correction and write brain maps
    if strcmp(mode, 'harvard-oxford')
        adj_p = zeros(2, length(stable_lr_rois));
        corrName = {'PearCorr', 'SpearCorr'};
        for mdl_idx = 1:2
            [h, crit_p, adj_ci_cvrg, adj_p(mdl_idx, :)]=fdr_bh(p_value(mdl_idx,:));
            
            p_mask = fill_mask_voxels_with_roi_values(lr_rois, p_value(mdl_idx, :));
            write_brain(p_mask, [], [], fullfile(output_dir, ...
                sprintf('pVal_%s_%s_stable%d_perm%d.nii', wordModel, corrName{mdl_idx}, NstableVox, Nperm)));
            p_corr_mask = fill_mask_voxels_with_roi_values(stable_lr_rois, adj_p(mdl_idx, :));
            write_brain(p_corr_mask, [], [], fullfile(output_dir, ...
                sprintf('pVal_fdr_%s_%s_stable%d_perm%d.nii', wordModel, corrName{mdl_idx}, NstableVox, Nperm)));

        end
        % save results
        save(sprintf(fullfile(output_dir, 'results_stable%d_perm%d.mat'), NstableVox, Nperm), 'p_value', 'adj_p', 'corrRDM');
    else 
        p_mask = fill_mask_voxels_with_roi_values({stable_lr_rois{1}}, 1);
            write_brain(p_mask, [], [], fullfile(output_dir, ...
                sprintf('wholebrain_map_%s_%s_stable%d_perm%d.nii', wordModel, 'SpearCorr', NstableVox, Nperm)));
        % save results
        save(sprintf(fullfile(output_dir, 'results_wholebrain_stable%d_perm%d.mat'), NstableVox, Nperm), 'p_value', 'corrRDM');
    end
    
                   
end