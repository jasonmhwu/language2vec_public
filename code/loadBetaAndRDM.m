function [beta, RDM, stimEntries] = loadBetaAndRDM(sessions, mode, wordModel, NstableVox)
% create ROI based on:
% all 96 ROIs
% specific ROIs
% specific ROIs + fusion
% can append word2vec representations at the end
% also outputs RDM
% sort words based on taxonomic then thematic categories


if nargin < 1
    sessions = subj_info();
end
if nargin < 2
    mode = 'harvard-oxford';
end
if nargin < 3
    wordModel = 'word2vec'; % 'none' means not appending word2vec
end
if nargin < 4
    NstableVox = 100;
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


[lr_rois, names] = load_rois(mode);


output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'compareRDM');
stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');
    
beta = cell(length(sessions), length(lr_rois));
RDM = cell(length(sessions), length(lr_rois));

for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    [all_betas_mat, runInfo, mask_mat] = ...
            retrieve_actualTrials_one_sub(session);
        
    output_dir = sprintf(output_dir_template, session.SubjectNum);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % standardize images
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
    if ~strcmp(mode, 'whole-brain')
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
        stimMatrix = sortrows(stimMatrix, [2, 3]);
        
        if strcmp(wordModel, 'word2vec')
            [~, stimMatrix] = load_wordModel_RDM(stimMatrix, wordModel);
            stimMatrix.beta = [stimMatrix.beta, stimMatrix.wordModel_vec];
            if EmptyfMRIFlag 
                stimMatrix.beta = stimMatrix.wordModel_vec;
            end
        end
        
        % output results
        beta{sess_idx, roi_idx} = stimMatrix.beta;
        RDM{sess_idx, roi_idx} = corr(stimMatrix.beta');
        
        
    end
                   
end
stimEntries = stimMatrix(:, [2, 3, 4]);