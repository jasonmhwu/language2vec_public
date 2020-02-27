function CreateGroupMask(sessions, mode, NstableVox, space, wordModel)
% can simple addition be a good model for specific areas?
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

if nargin < 1
    sessions = subj_info();
end
if nargin < 2  % 'harvard-oxford' or 'whole-brain'
    mode = 'whole-brain';
end
if nargin <3
    NstableVox = 100;
end
if nargin < 4 % 'beta' or 'similarity' or 'simDiff'
    space = 'beta';
end
if nargin < 5 % 'word2vec' or 'none'
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


if strcmp(mode, 'harvard-oxford')
    %rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
    %                                 'harvard_oxford_animal_space.nii');
    load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
    lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
elseif strcmp(mode, 'fedorenko-rois')
    load(fullfile(BaseDir(), 'fedorenko_rois', 'fedorenko_rois_realigned_all_rois.mat'));
    lr_rois = {roi_masks_struct.left};
end

output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'vectorAnalysis');
stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');
  
% initialize variables
Nsubj_intersect = 5;
EachSubjMask = cell(length(sessions), length(lr_rois));
EachSubj_maskmat = cell(1, length(sessions));
EachSubj_betas = cell(1, length(sessions));
GroupMask_union = cell(1, length(lr_rois));
GroupMask_intersect = cell(1, length(lr_rois));
session = sessions(1);
    [~, ~, mask_mat] = ...
            assemble_betas_one_sub(session.SubjectNum, 0, 0, 'SPM_analysis');
for roi = 1:length(lr_rois)
    GroupMask_union{roi} = zeros(size(mask_mat));
end

for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    [all_betas_mat, all_betas_info, mask_mat] = ...
            assemble_betas_one_sub(session.SubjectNum, 0, 0, 'SPM_analysis');
        EachSubj_maskmat{sess_idx} = mask_mat;
    runInfo = parseAndCombineOneSession(session.SubjectNum);
    
    % standardize images for each run
    runIdx = getColumn(runInfo, 'runNum');
    runIdx_uniq = unique(runIdx);
    
    for i = 1:length(runIdx_uniq)   
        tmp = all_betas_mat(find(runIdx == runIdx_uniq(i)), :);
        tmp = (tmp - mean(tmp)) ./ std(tmp);
        all_betas_mat(find(runIdx == runIdx_uniq(i)), :) = tmp;
    end
    EachSubj_betas{sess_idx} = all_betas_mat;
    % create stimulus matrix
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
            EachSubjMask{sess_idx, i} = lr_rois{i} & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
            GroupMask_union{i} = GroupMask_union{i} + double(EachSubjMask{sess_idx, i});
        end
    else % store whole brain as one ROI
        stabilityInRoi =sort(stableVoxels(logical(mask_mat)), 'descend');
        stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
        EachSubjMask{sess_idx, 1} = mask_mat & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
    end
end
% group subject masks together
for roi = 1:length(lr_rois)
    GroupMask_intersect{roi} = GroupMask_union{roi} >= Nsubj_intersect;
    GroupMask_union{roi} = GroupMask_union{roi} > 0;
end

% find overall 
rankResult = cell(1, length(lr_rois));
p_value = zeros(4, length(lr_rois));
for roi_idx = 1:length(lr_rois)
    EachSubj_roi_mask_indices = cell(length(sessions), length(lr_rois));
    for subj = 1:length(sessions)
        roi_mask_indices{subj, roi_idx} = rois_to_mask_indices( ...
            GroupMask_intersect(roi_idx), EachSubj_maskmat{subj}, true);
    end

    % average across runs to create beta matrix
    Stimulus_uniq = unique(Stimulus);
    stimMatrix = struct();
    for i = 1:length(Stimulus_uniq)
        trials = strcmp(Stimulus, Stimulus_uniq{i});
        stimMatrix(i).beta = mean(GroupRepre(trials, :));
        stimMatrix(i).wordType = wordType(find(trials, 1));
        stimMatrix(i).stimCategory = stimCategory(find(trials, 1));
        stimMatrix(i).Stimulus = Stimulus_uniq{i};
    end
    stimMatrix = struct2table(stimMatrix);

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
end
% save results
save(sprintf(fullfile(BaseDir(), 'cache', 'GroupMask_%s_%s_stable%d_intersect%d.mat'), ...
    wordModel, mode, NstableVox, Nsubj_intersect), 'GroupMask_union', 'GroupMask_intersect');

a=1;
