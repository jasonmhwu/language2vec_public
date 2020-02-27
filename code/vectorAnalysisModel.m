
function vectorAnalysisModel(sessions, mode, Nperm, NstableVox, space, wordModel)
% can simple addition be a good model for specific areas?
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

if nargin < 1
    sessions = subj_info();
end
if nargin < 2  % 'harvard-oxford' or 'whole-brain'
    mode = 'whole-brain';
end
if nargin < 3
    Nperm = 1000;
end
if nargin <4
    NstableVox = 100;
end
if nargin < 5 % 'beta' or 'similarity' or 'simDiff'
    space = 'beta';
end
if nargin < 6 % 'word2vec' or 'none'
    wordModel = 'none';
end

EmptyfMRIFlag = false;
if NstableVox == 0
    if strcmp(wordModel, 'none')
        error("you should always append word vectors when using zero fMRI voxels\n")';
    else
        NstableVox = 100;
        EmptyfMRIFlag = true;
    end
end

% load rois
[lr_rois, names] = load_rois(mode);

output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'vectorAnalysis');
stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');

for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    
    % specify and create directory to store results
    output_dir = sprintf(output_dir_template, session.SubjectNum);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    
    [all_betas_mat, runInfo, mask_mat] = ...
            retrieve_actualTrials_one_sub(session);
    
    wordType = getColumn(runInfo, 'wordType');
    Stimulus = getColumn(runInfo, 'Stimulus');
    stimCategory = getColumn(runInfo, 'stimCategory');
    
    % find the most stable voxels
    stable_fname = sprintf(stableVoxel_fname_template, session.SubjectNum);
    stableVoxels = spm_read_vols(spm_vol(stable_fname));
    if ~strcmp(mode, 'whole-brain')
        stable_lr_rois = cell(1, length(lr_rois));
        for i = 1:length(lr_rois)
            stabilityInRoi =sort(stableVoxels(lr_rois{i}), 'descend');
            stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
            stable_lr_rois{i} = lr_rois{i} & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
            % temporary: random voxels?
            %randVox = randsample(find(lr_rois{i}), min(NstableVox, length(stabilityInRoi)));
            %stable_lr_rois{i} = zeros(size(lr_rois{i}));
            %stable_lr_rois{i}(randVox) = 1;
        end
        if contains(mode, 'aggregate')
            tmp = zeros(size(lr_rois{1}));
            for i = 1:length(lr_rois)
                tmp = tmp + stable_lr_rois{i};
            end
            stable_lr_rois = cell(1);
            stable_lr_rois{1} = tmp;
        end
    else % store whole brain as one ROI
        stabilityInRoi =sort(stableVoxels(logical(mask_mat)), 'descend');
        stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
        stable_lr_rois{1} = mask_mat & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
    end
    
    roi_mask_indices = rois_to_mask_indices(stable_lr_rois, mask_mat);
    rankResult = cell(1, length(stable_lr_rois));
    p_value = zeros(4, length(stable_lr_rois));
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
        stimMatrix = sortrows(stimMatrix, [2, 3]);
        % append word2vec 
        if strcmp(wordModel, 'word2vec')
            [~, stimMatrix] = load_wordModel_RDM(stimMatrix, wordModel);
            stimMatrix.beta = [stimMatrix.beta, stimMatrix.wordModel_vec];
            if EmptyfMRIFlag 
                stimMatrix.beta = stimMatrix.wordModel_vec;
            end
        end
            
        categoryCombination = nchoosek(1:length(unique(stimCategory)), 2);
        typeCombination = nchoosek(unique(wordType), 2);

        % for one single pair of category * type combination, form a question,
        % and record ranks
        if strcmp(space, 'beta') || strcmp(space, 'similarity')
            rankResult{roi_idx} = rankCorrelation(stimMatrix, categoryCombination, typeCombination, space);
        else
            rankResult{roi_idx} = rankSimDiff(stimMatrix, categoryCombination, typeCombination);
        end
        % run permutation test
        % depending on what I permute, I can investigate whether
        % stimCategory or wordType is important for the model
        permResult = zeros(Nperm, length(rankResult{roi_idx}));
        parfor i = 1:Nperm
            tmpStim = stimMatrix;
            tmpStim.beta = tmpStim.beta(randperm(size(tmpStim, 1)), :);
            if strcmp(space, 'beta') || strcmp(space, 'similarity')
                result = rankCorrelation(tmpStim, categoryCombination, typeCombination, space);
            else
                result = rankSimDiff(tmpStim, categoryCombination, typeCombination);
            end
            tmpPermResult = zeros(1, length(result));
            for j = 1:length(result)
                tmpPermResult(j) = mean(result(j).meanRelativeRank(:));
            end
            permResult(i, :) = tmpPermResult;
        end
        
        % store permResults
        save(sprintf(fullfile(output_dir, 'permResult_%s_%s_stable%d_perm%d_subj%d.mat'), wordModel, mode, NstableVox, Nperm, sess_idx), 'permResult');
        
        realResult = zeros(1, length(rankResult{roi_idx}));
        for i = 1:length(rankResult{roi_idx})
            realResult(i) = mean(rankResult{roi_idx}(i).meanRelativeRank(:));
            p_value(i, roi_idx) = mean(realResult(i) < permResult(:, i));
            rankResult{roi_idx}(i).permPValue = p_value(i, roi_idx);
        end
        % report results
        if strcmp(mode, 'harvard-oxford')
            load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
        
            if roi_idx <= 48
                fprintf('left %s has %d voxels.\n', roi_masks_struct(roi_idx).name, size(beta_masked, 2));
            else
                fprintf('right %s has %d voxels.\n', roi_masks_struct(roi_idx-48).name, size(beta_masked, 2));
            end
            fprintf('rankAllWords model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(1), 0.5, p_value(1, roi_idx));
            fprintf('rankWithinType model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(2), 0.5, p_value(2, roi_idx));
            fprintf('taxonomic rank model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(3), 0.5, p_value(3, roi_idx));
             fprintf('thematic rank model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(4), 0.5, p_value(4, roi_idx));
        elseif strcmp(mode, 'fedorenko-rois')
             fprintf('left %s has %d voxels.\n', names{roi_idx}, size(beta_masked, 2));
             fprintf('rankAllWords model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(1), 0.5, p_value(1, roi_idx));
            fprintf('rankWithinType model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(2), 0.5, p_value(2, roi_idx));
            fprintf('taxonomic rank model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(3), 0.5, p_value(3, roi_idx));
             fprintf('thematic rank model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(4), 0.5, p_value(4, roi_idx));
        elseif strcmp(mode, 'whole-brain')
            fprintf('most stable whole brain has %d voxels.\n', size(beta_masked, 2));
            fprintf('rankAllWords model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(1), 0.5, p_value(1));
            fprintf('rankWithinType model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(2), 0.5, p_value(2));
            fprintf('rankTypeAgainstType model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(3), 0.5, p_value(3));
            fprintf('thematic rank model have an meanRelative rank of %f, where chance is %f, permutation test p-value is %f\n', ...
                realResult(4), 0.5, p_value(4));
        else
        end
    end % end roi_idx loop
    
    
    
    if EmptyfMRIFlag % for storing purposes
        NstableVox = 0;
    end
    % fdr correction and write brain maps
    if ~strcmp(mode, 'whole-brain') & ~contains(mode, 'aggregate')
        meanRank = zeros(4, length(stable_lr_rois));
        for mdl_idx = 1:4
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_value(mdl_idx,:));
            for i = 1:length(stable_lr_rois)
                rankResult{i}(mdl_idx).p_corr = adj_p(i);
            end
            p_mask = fill_mask_voxels_with_roi_values(lr_rois, -log(p_value(mdl_idx, :)));
            write_brain(p_mask, [], [], fullfile(output_dir, ...
                sprintf('pVal_%s_stable%d_perm%d.nii', rankResult{1}(mdl_idx).modelName, NstableVox, Nperm)));
            p_corr_mask = fill_mask_voxels_with_roi_values(lr_rois, -log(adj_p));
            write_brain(p_corr_mask, [], [], fullfile(output_dir, ...
                sprintf('pVal_fdr_%s_stable%d_perm%d.nii', rankResult{1}(mdl_idx).modelName, NstableVox, Nperm)));
            % write mean rank map
            for i = 1:length(stable_lr_rois)
                meanRank(mdl_idx, i) = mean(rankResult{i}(mdl_idx).meanRelativeRank(:));
                rankResult{i}(mdl_idx).MeanRank = mean(rankResult{i}(mdl_idx).meanRelativeRank(:));
            end
            rank_mask = fill_mask_voxels_with_roi_values(lr_rois, meanRank(mdl_idx, :));
            write_brain(rank_mask, [], [], fullfile(output_dir, ...
                sprintf('rank_%s_stable%d_perm%d.nii', rankResult{1}(mdl_idx).modelName, NstableVox, Nperm)));
            
        end
    end   
    
    % save results
        save(sprintf(fullfile(output_dir, 'results_%s_%s_stable%d_perm%d.mat'), wordModel, mode, NstableVox, Nperm), 'rankResult');
     
        
     if EmptyfMRIFlag % for later processing purposes
        NstableVox = 100;
    end
end
end


function [rankResult, hat_true, hat_wrong] = rankCorrelation(stimMatrix, categoryCombination, typeCombination, space)
% rank the correlation between model word and all other words in the
% stimulus list
    rankResult = struct();
    rankResult(1).modelName = 'predict-identity';
    rankResult(2).modelName = 'rankWithinType';
    rankResult(3).modelName = 'predict-category';
    rankResult(4).modelName = 'predict-theme';
    %rankResult(1).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    %rankResult(2).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    %rankResult(3).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    %rankResult(4).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    
    % temporary:
    %rankResult(2).modelName = 'rankAllWords';
    %rankResult(3).modelName = 'rankFourRelevantWords';
    %rankResult(4).modelName = 'rankSameTheme';
    %rankResult(5).modelName = 'rankSameCategory';
    %rankResult(5).meanRelativeRank = nan(size(categoryCombination, 1), size(typeCombination, 1));
    
    hat_true = [];
    hat_wrong = [];
    
    
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
           % option to use beta space or similarity space
           if strcmp(space, 'beta')
                Betas = stimMatrix.beta;
           else
                Betas = corr(stimMatrix.beta');
                wordDim = 1:size(stimMatrix, 1);
                wordDim = wordDim(~ismember(wordDim, words_idx));
                Betas = Betas(:, wordDim);
           end
           combinedBeta = [Betas(words_idx(2), :) + Betas(words_idx(3), :) - Betas(words_idx(4), :); ...
                                       Betas(words_idx(1), :) + Betas(words_idx(4), :) - Betas(words_idx(3), :); ...
                                       Betas(words_idx(1), :) + Betas(words_idx(4), :) - Betas(words_idx(2), :); ...
                                       Betas(words_idx(2), :) + Betas(words_idx(3), :) - Betas(words_idx(1), :)];
           
           % for each question, which word to compute correlation and rank
           % with
           correlations = corr(Betas', combinedBeta');
           
           % predict-identity model
           % theoretical chance: 0.5
           wordsToCompare = ~ismember(1:size(stimMatrix, 1), words_idx);
           rank1 = zeros(1, 4); % 4 words to test in each analogy
           for i = 1:4
               rank1(i) = sum(correlations(words_idx(i), i) > correlations(wordsToCompare, i)) / ...
                   sum(wordsToCompare);
           end
           rankResult(1).meanRelativeRank(cat_idx, type_idx) = mean(rank1);
           
           % temporary model 2: compare true word to all words
           %{
           rank2 = zeros(1, 4);
           for i = 1:4
               rank2(i) = sum(correlations(words_idx(i), i) > correlations(:, i)) / ...
                   (size(correlations, 1) - 1);
           end
           rankResult(2).meanRelativeRank(cat_idx, type_idx) = mean(rank2);
           
           % temporary: compare true word to the four words
           rank3 = zeros(1, 4);
           for i = 1:4
               rank3(i) = sum(correlations(words_idx(i), i) > correlations(words_idx, i)) / ...
                   3;
           end
           rankResult(3).meanRelativeRank(cat_idx, type_idx) = mean(rank3);
           
           % temporary: compare true word to all words in the same theme
           rank4 = zeros(1, 4);
           for i = 1:4
               wordsToCompare = find(stimMatrix.stimCategory == stimMatrix.stimCategory(words_idx(i)));
               rank4(i) = sum(correlations(words_idx(i), i) > correlations(wordsToCompare, i)) / ...
                   2;
           end
           
           rankResult(4).meanRelativeRank(cat_idx, type_idx) = mean(rank4);
            %}

           % rankWithinType model
           % predict theme
           % theoretical chance: 0.5
           wordsToCompare = ~ismember(1:size(stimMatrix, 1), words_idx);
           
           rank2 = zeros(1, 4);
           for i = 1:4
               sameType = strcmp(stimMatrix.wordType, stimMatrix.wordType(words_idx(i)))';
               rank2(i) = sum(correlations(words_idx(i), i) > correlations(wordsToCompare & sameType, i)) / ...
                   (sum(wordsToCompare & sameType));
           end
           rankResult(2).meanRelativeRank(cat_idx, type_idx) = mean(rank2);
           
           
           % rankTypeAgainstType model
           % taxonomic rank model
           % theoretical chance: 0.5
           % are candidate words of true word type higher ranked than other
           % word types?
           rank3 = zeros(1, 4);
           % find the 4 words
           TypeLastCategory = {'tool'; 'person'; 'building'};
           words_idx_last_category = [find(stimMatrix.stimCategory == categoryCombination(cat_idx, 1) & ...
                         strcmp(stimMatrix.wordType, TypeLastCategory{type_idx})), ...
                         find(stimMatrix.stimCategory == categoryCombination(cat_idx, 2) & ...
                         strcmp(stimMatrix.wordType, TypeLastCategory{type_idx}))];
           wordsToCompare = ~ismember(1:size(stimMatrix, 1), ...
               [words_idx, words_idx_last_category]);
           
           for i = 1:4
               
               sameType = strcmp(stimMatrix.wordType, ...
                   stimMatrix.wordType(words_idx(i)))';
               sameType = find(sameType & wordsToCompare);
               otherTypeName = TypeLastCategory( ...
                   ~contains(TypeLastCategory, stimMatrix.wordType(words_idx(i))));
               diffType = {find(strcmp(stimMatrix.wordType, otherTypeName{1})' & wordsToCompare); ...
                                  find(strcmp(stimMatrix.wordType, otherTypeName{2})' & wordsToCompare)};
               
               sameTypeTemplate = mean(stimMatrix.beta(sameType, :), 1);
               diffTypeTemplate = [mean(stimMatrix.beta(diffType{1}, :), 1); ...
                   mean(stimMatrix.beta(diffType{2}, :), 1)];
               
               rank3(i) = sum(corr(combinedBeta(i, :)', sameTypeTemplate') > ...
                   corr(combinedBeta(i, :)', diffTypeTemplate')) / 2;
                            
               % obsolete versions
               %{
               sameType = strcmp(stimMatrix.wordType, stimMatrix.wordType(words_idx(i)))';
               diffType = ~sameType;
               sameType = find(sameType & wordsToCompare);
               rank_individual = zeros(1, length(sameType));
               for j = 1:length(sameType)
                   % temporary: try to only compare words in the same theme
                   sameTheme = (stimMatrix.stimCategory == stimMatrix.stimCategory(sameType(j)))';
                   rank_individual(j) = sum(correlations(sameType(j), i) > ...
                        correlations(sameTheme & diffType, i)) / ...
                        (sum(sameTheme & diffType));
                   
                    % original: compare all 26 words from the different
                    % type
                    %rank_individual(j) = sum(correlations(sameType(j), i) > ...
                     %   correlations(wordsToCompare & diffType, i)) / ...
                      %  (sum(wordsToCompare & diffType));
               end
               
               rank3(i) = mean(rank_individual);
               %}
           end
           rankResult(3).meanRelativeRank(cat_idx, type_idx) = mean(rank3);
           
           % thematic rank model
           % theoretical chance: 0.5
           % are candidate words of true thematic category higher ranked than 
           % others?
           rank4 = zeros(1, 4);
           % find the 4 words
           TypeLastCategory = {'tool'; 'person'; 'building'};
           words_idx_last_category = [find(stimMatrix.stimCategory == categoryCombination(cat_idx, 1) & ...
                         strcmp(stimMatrix.wordType, TypeLastCategory{type_idx})), ...
                         find(stimMatrix.stimCategory == categoryCombination(cat_idx, 2) & ...
                         strcmp(stimMatrix.wordType, TypeLastCategory{type_idx}))];
           wordsToCompare = strcmp(stimMatrix.wordType, TypeLastCategory{type_idx});
           wordsToCompare(words_idx_last_category) = false;
           
           
           for i = 1:2
               rank4(i)= sum(correlations(words_idx_last_category(1), i) > ...
                   correlations(wordsToCompare, i)) / ...
                        (sum(wordsToCompare));
           end
           
           for i = 3:4
               rank4(i)= sum(correlations(words_idx_last_category(2), i) > ...
                   correlations(wordsToCompare, i)) / ...
                        (sum(wordsToCompare));
           end
           rankResult(4).meanRelativeRank(cat_idx, type_idx) = mean(rank4);
           
           
        end   
    end
end

function rankResult = rankSimDiff(stimMatrix, categoryCombination, typeCombination)
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
           % option to use beta space or similarity space
           
            Betas = corr(stimMatrix.beta'); % actually similarity
           % for each question, which word to compute correlation and rank
           % with
           wordsToCompare = [~ismember(1:size(stimMatrix, 1), words_idx([2, 3, 4])); ...
                             ~ismember(1:size(stimMatrix, 1), words_idx([1, 3, 4])); ...
                             ~ismember(1:size(stimMatrix, 1), words_idx([1, 2, 4])); ...
                             ~ismember(1:size(stimMatrix, 1), words_idx([1, 2, 3]))];
           SimAbsDiff = [abs(Betas(2, :) - Betas(4, 3)) + abs(Betas(3, :) - Betas(2, 4));
               abs(Betas(1, :) - Betas(4, 3)) + abs(Betas(4, :) - Betas(1, 3));
               abs(Betas(1, :) - Betas(4, 2)) + abs(Betas(4, :) - Betas(1, 2));
               abs(Betas(3, :) - Betas(1, 2)) + abs(Betas(2, :) - Betas(1, 3))];
               
           
           % rankAllWords model
           % theoretical chance: 0.5
           rank1 = zeros(1, 4);
           for i = 1:4
               rank1(i) = sum(SimAbsDiff(i, words_idx(i)) < SimAbsDiff(i, wordsToCompare(i, :))) / ...
                   sum(wordsToCompare(i, :));
           end
           rankResult(1).meanRelativeRank(cat_idx, type_idx) = mean(rank1);
           
           % rankWithinType model
           % theoretical chance: 0.5
           rank2 = zeros(1, 4);
           for i = 1:4
               sameType = strcmp(stimMatrix.wordType, stimMatrix.wordType(words_idx(i)))';
               rank2(i) = sum(SimAbsDiff(i, words_idx(i)) < SimAbsDiff(i, wordsToCompare(i, :) & sameType)) / ...
                   sum(wordsToCompare(i, :) & sameType);
           end
           rankResult(2).meanRelativeRank(cat_idx, type_idx) = mean(rank2);
        end   
    end
end