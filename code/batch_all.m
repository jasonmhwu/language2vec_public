%%% Run the whole preprocessing pipeline on everything
subjects = subj_info();
%subjects = subjects(1:(end-1));
new_subj = subjects(end);
batch_realign(new_subj);
batch_sliceTimeCorrect(new_subj);
batch_coregister(new_subj);
batch_segment(new_subj);
batch_normalize_write(new_subj);
batch_model_spec_and_est(subjects, BaseDir(), 'SPM_analysis_motionReg');
batch_contrasts(new_subj);

% include motion regressors in this version

% second branch - smoothed data, only used for glm contrast!!
batch_smooth(new_subj);
batch_model_spec_and_est(new_subj, BaseDir(), ...
    'smoothed_SPM_analysis/', 'swar');
batch_contrasts(new_subj, 'smoothed_SPM_analysis/');

% stable voxels analysis
parpool(8);
stableVoxels(subjects);
stableVoxels_location(subjects, 1600);
% three-way decoding 
% normal rad=3 searchlight
do_multignb_wordtype(subj_info(),  'multisvm_wordtype_zscore');
% perform on 800 whole-brain voxels
do_multignb_main(subjects,  'multisvm_theme_whole-brain', 800);
do_multignb_main(subjects,  'multisvm_category_whole-brain', 800);

multisvm_group_level('whole-brain', 'theme', 800);

% perform on 96 Harvard-Oxford Atlases
% do group-level analysis to see if I can find significant ROIs
% hopefully use these ROIs to justify further vecAnalysis
voxSizeList = [50, 100, 200, 400, 800];
do_multignb_main(subjects,  'multisvm_theme_harvard-oxford', voxSizeList);
do_multignb_main(subjects,  'multisvm_category_harvard-oxford', voxSizeList);

multisvm_group_level('harvard-oxford', 'category', voxSizeList(1));

%% vector analysis - harvard oxford
%voxSizeList =  [0, 100, 200, 400, 800, 1600];
voxSizeList = [50];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'harvard-oxford', 2, voxSizeList(voxSize_idx), 'beta', 'word2vec');
    vectorAnalysis_group(subjects, 'harvard-oxford', 2, voxSizeList(voxSize_idx), 'beta', 'word2vec');
end
%voxSizeList =  [100, 200, 400, 800, 1600];
voxSizeList = [50];
for voxSize_idx = 1:length(voxSizeList)
    %vectorAnalysisModel(subjects, 'harvard-oxford', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'harvard-oxford', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end

%% vector analysis - fedorenko_rois
voxSizeList =  [0, 50, 100, 200, 400, 800];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'fedorenko-rois', 2, voxSizeList(voxSize_idx), 'beta', 'word2vec');
    vectorAnalysis_group(subjects, 'fedorenko-rois', 2, voxSizeList(voxSize_idx), 'beta', 'word2vec');
end
voxSizeList =  [400];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'fedorenko-rois', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'fedorenko-rois', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end
vectorAnalysis_group_plotAcc_ROI('fedorenko-rois', 2, voxSizeList, 'beta');

%% vector analysis - MNI four lobes
voxSizeList = [50, 100, 200, 400];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'MNI-four-lobes', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'MNI-four-lobes', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end
vectorAnalysis_group_plotAcc_ROI('MNI-four-lobes', 2, voxSizeList, 'beta');
%% vector analysis - MNI four lobes - aggregated
voxSizeList = [25, 50, 100, 200, 400];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'MNI-four-lobes-aggregate', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'MNI-four-lobes-aggregate', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end
vectorAnalysis_group_plotAcc_ROI('MNI-four-lobes-aggregate', 2, voxSizeList, 'beta');

%% vector analysis - whole brain
voxSizeList =  [0, 100, 200, 400, 800, 1600];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'whole-brain', 2, voxSizeList(voxSize_idx), 'beta', 'word2vec');
    vectorAnalysis_group(subjects, 'whole-brain', 2, voxSizeList(voxSize_idx), 'beta', 'word2vec');
end
%voxSizeList = 100:100:1600;
voxSizeList = [800];
for voxSize_idx = 1:length(voxSizeList)
    %vectorAnalysisModel(subj_info(), 'whole-brain', 10000, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysisPermTest(subj_info(), 'whole-brain', 10000, voxSizeList(voxSize_idx), 'beta', 'none');
 
    %vectorAnalysis_group(subjects, 'whole-brain', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end
vectorAnalysis_group_plotAcc('whole-brain', 2, voxSizeList, 'beta');

%% vector analysis - union of fedorenko vs. other part of the whole brain
voxSizeList =  [800];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'union-language-network', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'union-language-network', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end

%% vector analysis - general semantic network by Jeffrey Binder meta-analysis
voxSizeList =  100;
for voxSize_idx = 1:length(voxSizeList)
    %vectorAnalysisModel(subjects, 'general-semantic-network', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'general-semantic-network', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end

%% vector analysis - general semantic network by Jeffrey Binder meta-analysis
voxSizeList =  [10, 20, 40, 80, 160];
for voxSize_idx = 2:length(voxSizeList)
    vectorAnalysisModel(subjects, 'general-semantic-network-aggregate', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'general-semantic-network-aggregate', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end
%% vector analysis - cortical and subcortical Harvard-Oxford
voxSizeList =  [800];
for voxSize_idx = 1:length(voxSizeList)
    vectorAnalysisModel(subjects, 'harvard-oxford-cortical-subcortical', 2, voxSizeList(voxSize_idx), 'beta', 'none');
    vectorAnalysis_group(subjects, 'harvard-oxford-cortical-subcortical', 2, voxSizeList(voxSize_idx), 'beta', 'none');
end

%% RDM Analysis - whole brain
for voxSize = [200, 400, 800]
    RDMAnalysis(subjects, 'whole-brain', 'word2vec', voxSize);
    RDMAnalysis(subjects, 'whole-brain', 'none', voxSize);
end
RDMAnalysis(subjects, 'whole-brain', 'word2vec', 0); % pure word2vec vectors, no fMRI
RDMAnalysis_group_plotAcc('whole-brain', [200, 400, 800]);

%% RDM Analysis - harvard oxford
for voxSize = [200, 400, 800]
    RDMAnalysis(subjects, 'harvard-oxford', 'word2vec', voxSize);
    RDMAnalysis(subjects, 'harvard-oxford', 'none', voxSize);
end

%% RDM Analysis - fedorenko-rois
for voxSize = [50, 100, 200, 400, 800]
    RDMAnalysis(subjects, 'fedorenko-rois', 'word2vec', voxSize);
    RDMAnalysis(subjects, 'fedorenko-rois', 'none', voxSize);
end
RDMAnalysis(subjects, 'fedorenko-rois', 'word2vec', 0);
RDMAnalysis_group_plotAcc('fedorenko-rois', [50, 100, 200, 400, 800]);

%% RDM Difference Analysis
for voxSize = [200]
    RDMDifferenceAnalysis(subjects, 'whole-brain', 'none', voxSize);
end

%% Create Group Mask
CreateGroupMask(subjects, 'fedorenko-rois', 200, 'beta', 'word2vec');

%% run stuff with group mask

%% estimate noise ceiling
estimateNoiseCeiling(subjects, 'harvard-oxford', 200);

% MRI quality check
checkHeadMotion;

%% obselete analyses
% vector analysis on similarity vector
%vectorAnalysisModel(subjects, 'whole-brain', 1000, 800, 'similarity');

% find the closest word with similar similarity
%vectorAnalysisModel(subjects, 'whole-brain', 1000, 800, 'simDiff');

% compare word2vec and roi RDM 
compareRDM_word2vec_image(subjects, 'harvard-oxford', 'word2vec', 200, 1000);
compareRDM_word2vec_image(new_subj, 'whole-brain', 'word2vec', 50, 500);

% intermediate function: load betas and RDM for further t-test or RSA
% can append word2vec representation here
% [beta, RDM] = loadBetaAndRDM(subjects, 'whole-brain', 'word2vec', 200);