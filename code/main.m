% main script to replicate paper results

% load information about 13 subjects
subjects = subj_info();

% assign each voxel a stability score
% stability score is computed as the summation of cross-run correlations
% results will be stored in subjects/Subject%03d/stableVoxels
stableVoxels(subjects);

%% replicate result 1
% take most stable 800 voxels out of whole-brain, and apply 3 ranking
% metrics
voxSize = 800;
% performs the 3 ranking metrics
vectorAnalysisModel(subjects, 'whole-brain', 2, voxSize, 'beta', 'none');
% aggregate group-level results and do t-test
vectorAnalysis_group(subjects, 'whole-brain', 2, voxSize, 'beta', 'none');

%% replicate result 2
% take most stable 100 voxels from each of 8 semantic ROIs, and apply 3 ranking
% metrics
voxSize = 100;
vectorAnalysisModel(subjects, 'general-semantic-network', 2, voxSize, 'beta', 'none');
vectorAnalysis_group(subjects, 'general-semantic-network', 2, voxSize, 'beta', 'none');
