function estimateNoiseCeiling(sessions, mode, NstableVox)
% sum up 
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

if nargin < 1
    sessions = subj_info();
end
if nargin < 2
    mode = 'harvard-oxford';
end
if nargin < 3
    NstableVox = 100;
end


if strcmp(mode, 'harvard-oxford')
    %rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
    %                                 'harvard_oxford_animal_space.nii');
    load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
    lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
end

output_dir_template = fullfile(BaseDir(), 'noiseCeiling');
stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');
   
subject_roi_RDM = cell(1, length(sessions));

for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    [all_betas_mat, all_betas_info, mask_mat] = ...
            assemble_betas_one_sub(session.SubjectNum, 0, 0, 'SPM_analysis');
        
    output_dir = sprintf(output_dir_template, session.SubjectNum);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % standardize images, one for each run
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
    roi_RDM = cell(1, length(stable_lr_rois));
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
        
        % calculate real fMRI RDM
        roi_RDM{roi_idx} = 1 - corr(stimMatrix.beta');
    end
    
    subject_roi_RDM{sess_idx} = roi_RDM;
                   
end

% estimate noise ceiling
lowerNoiseCeil = zeros(length(sessions), length(lr_rois));
upperNoiseCeil = zeros(length(sessions), length(lr_rois));
roi_idx = 31;
for sess_idx = 1:length(sessions)
    %for roi_idx = 1:length(lr_rois)
        
        sumRDM = zeros(size(subject_roi_RDM{1}{1}));
        for i = 1:length(sessions)
            if i ~= sess_idx
                sumRDM = sumRDM + subject_roi_RDM{i}{roi_idx};
            end
        end
        % lower noise ceiling
        lowerNoiseCeil(sess_idx, roi_idx) = corr(sumRDM(triu(true(size(sumRDM)),1)), ...
                                   subject_roi_RDM{sess_idx}{roi_idx}(triu(true(size(sumRDM)),1)), ...
                                   'type', 'Spearman');
        sumRDM = sumRDM + subject_roi_RDM{sess_idx}{roi_idx};
        % higher noise ceiling
        upperNoiseCeil(sess_idx, roi_idx) = corr(sumRDM(triu(true(size(sumRDM)),1)), ...
                                   subject_roi_RDM{sess_idx}{roi_idx}(triu(true(size(sumRDM)),1)), ...
                                   'type', 'Spearman');
   % end
end
imagesc(sumRDM);
% write brain image
mask = fill_mask_voxels_with_roi_values(lr_rois, mean(lowerNoiseCeil));
write_brain(mask, [], [], fullfile(output_dir, ...
    sprintf('lowerNoiseCeil_stable%d.nii', NstableVox)));
mask = fill_mask_voxels_with_roi_values(lr_rois, mean(upperNoiseCeil));
write_brain(mask, [], [], fullfile(output_dir, ...
    sprintf('upperNoiseCeil_stable%d.nii', NstableVox)));
a=1;
