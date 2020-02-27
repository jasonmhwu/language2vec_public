function similarityVectorAnalysisModel(sessions, mode, Nperm, NstableVox)
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
    Nperm = 1000;
end
if nargin <4
    NstableVox = 100;
end

if strcmp(mode, 'harvard-oxford')
    %rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
    %                                 'harvard_oxford_animal_space.nii');
    load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
    lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
end

output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'vectorAnalysis');
stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');
    
for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    [all_betas_mat, all_betas_info, mask_mat] = ...
            assemble_betas_one_sub(session.SubjectNum, 0, 0, 'SPM_analysis');
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
    rankResult = cell(1, length(stable_lr_rois));
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

        categoryCombination = nchoosek(1:length(unique(stimCategory)), 2);
        typeCombination = nchoosek(unique(wordType), 2);

        % for one single pair of category * type combination, form a question,
        % and record ranks
        rankResult{roi_idx} = rankCorrelation(stimMatrix, categoryCombination, typeCombination);
        
        % run permutation test
        % depending on what I permute, I can investigate whether
        % stimCategory or wordType is important for the model
        permResult = zeros(length(rankResult{roi_idx}), Nperm);
        parfor i = 1:Nperm
            tmpStim = stimMatrix;
            tmpStim.beta = tmpStim.beta(randperm(size(tmpStim, 1)), :);
            
            result = rankCorrelation(tmpStim, categoryCombination, typeCombination);
            permResult(:, i) = [mean(result(1).meanRelativeRank(:)), mean(result(2).meanRelativeRank(:))];
        end
        realResult = [mean(rankResult{roi_idx}(1).meanRelativeRank(:)), mean(rankResult{roi_idx}(2).meanRelativeRank(:))];
        p_value(:, roi_idx) = [mean(realResult(1) < permResult(1, :)), mean(realResult(2) < permResult(2, :))];
        rankResult{roi_idx}(1).permPValue = p_value(1, roi_idx);
        rankResult{roi_idx}(2).permPValue = p_value(2, roi_idx);
        % report results
        if strcmp(mode, 'harvard-oxford')
            if roi_idx <= 48
                fprintf('left %s has %d voxels.\n', roi_masks_struct(roi_idx).name, size(beta_masked, 2));
            else
                fprintf('right %s has %d voxels.\n', roi_masks_struct(roi_idx-48).name, size(beta_masked, 2));
            end
            fprintf('rankAllWords model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(1), mean(0:53)/54, p_value(1, roi_idx));
            fprintf('rankWithinType model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(2), mean(0:17)/18, p_value(2, roi_idx));
        else
            fprintf('most stable whole brain has %d voxels.\n', size(beta_masked, 2));
            fprintf('rankAllWords model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(1), mean(0:53)/54, p_value(1));
            fprintf('rankWithinType model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(2), mean(0:17)/18, p_value(2));
        end
    end
    
    output_dir = sprintf(output_dir_template, session.SubjectNum);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % fdr correction and write brain maps
    if strcmp(mode, 'harvard-oxford')
        for mdl_idx = 1:2
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_value(mdl_idx,:));
            for i = 1:length(stable_lr_rois)
                rankResult{roi_idx}(mdl_idx).p_corr = adj_p(i);
            end
            p_mask = fill_mask_voxels_with_roi_values(lr_rois, p_value(mdl_idx, :));
            write_brain(p_mask, [], [], fullfile(output_dir, ...
                sprintf('pVal_%s_stable%d_perm%d.nii', rankResult{1}(mdl_idx).modelName, NstableVox, Nperm)));
            p_corr_mask = fill_mask_voxels_with_roi_values(lr_rois, adj_p);
            write_brain(p_corr_mask, [], [], fullfile(output_dir, ...
                sprintf('pVal_fdr_%s_stable%d_perm%d.nii', rankResult{1}(mdl_idx).modelName, NstableVox, Nperm)));

        end
        % save results
        save(sprintf(fullfile(output_dir, 'results_stable%d_perm%d.mat'), NstableVox, Nperm), 'rankResult');
    else 
        % save results
        save(sprintf(fullfile(output_dir, 'results_wholebrain_stable%d_perm%d.mat'), NstableVox, Nperm), 'rankResult');
    end
    
                   
end

function rankResult = rankCorrelation(stimMatrix, categoryCombination, typeCombination)
% rank the correlation between model word and all other words in the
% stimulus list
    rankResult = struct();
    rankResult(1).modelName = 'rankAllWords';
    rankResult(2).modelName = 'rankWithinType';
    rankResult(1).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    rankResult(2).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    for cat_idx = 1:size(categoryCombination, 1)
        for type_idx = 1:size(typeCombination, 1)
            % find the 4 words
            words_idx = [find(stimMatrix.stimCategory == categoryCombination(cat_idx, 1) & ...
                         strcmp(stimMatrix.wordType, typeCombination{type_idx, 1})), ...
                         find(stimMatrix.stimCategory == categoryCombination(cat_idx, 1) & ...
                         strcmp(stimMatrix.wordType, typeCombination{type_idx, 2})), ...
                         find(stimMatrix.stimCategory == categoryCombination(cat_idx, 2) & ...
                         strcmp(stimMatrix.wordType, typeCombination{type_idx, 1})), ...
                         find(stimMatrix.stimCategory == categoryCombination(cat_idx, 2) & ...
                         strcmp(stimMatrix.wordType, typeCombination{type_idx, 2}))];
                     
           % they form four questions, now compute the synthesized word
           Betas = stimMatrix.beta;
           combinedBeta = [Betas(words_idx(2), :) + Betas(words_idx(3), :) - Betas(words_idx(4), :); ...
                           Betas(words_idx(1), :) - Betas(words_idx(3), :) + Betas(words_idx(4), :); ...
                           Betas(words_idx(1), :) - Betas(words_idx(2), :) + Betas(words_idx(4), :); ...
                          -Betas(words_idx(1), :) + Betas(words_idx(2), :) + Betas(words_idx(3), :)];
           % for each question, which word to compute correlation and rank
           % with
           wordsToCompare = [~ismember(1:size(stimMatrix, 1), words_idx([2, 3, 4])); ...
                             ~ismember(1:size(stimMatrix, 1), words_idx([1, 3, 4])); ...
                             ~ismember(1:size(stimMatrix, 1), words_idx([1, 2, 4])); ...
                             ~ismember(1:size(stimMatrix, 1), words_idx([1, 2, 3]))];
           
           
           correlations = corr(Betas', combinedBeta');
           % rankAllWords model
           % theoretical chance: 0.4907
           rank1 = zeros(1, 4);
           for i = 1:4
               rank1(i) = sum(correlations(words_idx(i), i) > correlations(wordsToCompare(i, :), i)) / ...
                   sum(wordsToCompare(i, :));
           end
           rankResult(1).meanRelativeRank(cat_idx, type_idx) = mean(rank1);
           
           % rankWithinType model
           % theoretical chance: 0.4722
           rank2 = zeros(1, 4);
           for i = 1:4
               sameType = strcmp(stimMatrix.wordType, stimMatrix.wordType(words_idx(i)))';
               rank2(i) = sum(correlations(words_idx(i), i) > correlations(wordsToCompare(i, :) & sameType, i)) / ...
                   sum(wordsToCompare(i, :) & sameType);
           end
           rankResult(2).meanRelativeRank(cat_idx, type_idx) = mean(rank2);
        end   
    end
