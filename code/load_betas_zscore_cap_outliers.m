% load_betas_zscore_cap_outliers
%
% loads single-trial beta images, extracts the ones for actual trials, z scores
% them (within each voxel), and caps outliers (sets them to the maximum allowed
% signed z value).
%
% Input: 
%   subject: subject info struct (passed to assemble_betas_one_subject)
%   [z_score_cutoff=4]: maximum allowed (absolute) z-score (optional)
%   [trial_exclude_outlier_proportion=Inf]: maximum allowed proportion of
%     outliers for a given trial before it's excluded (optional, default to NO
%     exclusions of trials)
%
% Output: 
%   beta_z: trial x voxel matrix of z-scored (and outlier capped) beta values.
%   beta_info: the nifti headers corresponding to the rows of beta_z.
%   mask_mat: the mask volume to map voxel indices back into xyz coordinates
%   excluded: a logical matrix indicating which of the beta_z entries were
%     identified as outliers
%   n_outliers_by_trial: number of outliers identified in each trial (use to
%     detect weird trials)
%   n_outliers_by_voxel: number of outliers identified in each voxel (use to
%     detect voxels with highly non-normal distributions, etc.)
%   trials_included: logical vector indicating which trials are included in the
%     output (beta_z and beta_info).  If outlier proportion cutoff is not
%     specified, this will be all trials present in the input.
function [beta_z, beta_info, mask_mat, outliers, n_outliers_by_trial, n_outliers_by_voxel, trials_include_mask] = ...
    load_betas_zscore_cap_outliers(subject, z_score_cutoff, trial_exclude_outlier_proportion)

persistent last_beta_z last_beta_info last_mask_mat last_outliers ...
    last_n_outliers_by_trial last_n_outliers_by_voxel last_trials_include_mask ...
    last_subject last_z_score_cutoff last_trial_exclude_outlier_proportion  

import functools.*;

% get rid of betas where abs z score is above z_score_cutoff (default 4)
if nargin < 2
    z_score_cutoff = 4;
end

if nargin < 3
    trial_exclude_outlier_proportion = Inf;
else
    if trial_exclude_outlier_proportion > 1 || trial_exclude_outlier_proportion < 0
        error('Proportion of outlier voxels cutoff is outside the [0, 1] range: %d', ...
            trial_exclude_outlier_proportion);
    end
end

%% check if persistent variables are all non-empty and match input
persistents_empty = map(@isempty, ...
                        {last_beta_z, last_beta_info, last_mask_mat, last_outliers, ...
                         last_n_outliers_by_trial, last_n_outliers_by_voxel, ...
                         last_trials_include_mask, last_subject, last_z_score_cutoff, ...
                         last_trial_exclude_outlier_proportion});
if ~any([persistents_empty{:}]) && ...
        last_subject.SubjectNum == subject.SubjectNum && ...
        strcmp(last_subject.PrePost, subject.PrePost) && ...
        last_z_score_cutoff == z_score_cutoff && ...
        last_trial_exclude_outlier_proportion == trial_exclude_outlier_proportion
    beta_z = last_beta_z;
    beta_info = last_beta_info;
    mask_mat = last_mask_mat;
    outliers = last_outliers;
    n_outliers_by_trial = last_n_outliers_by_trial;
    n_outliers_by_voxel = last_n_outliers_by_voxel;
    trials_include_mask = last_trials_include_mask;
    return
end

%% can't use cached values, do it live

[beta_mat, beta_info, mask_mat] = assemble_betas_one_sub(subject.SubjectNum, 1, 1);

% ignore constant predictors for purposes of calculating outliers.
non_baseline_rows = cell2mat(map(@(d) isempty(strfind(d, 'constant')), {beta_info.descrip}));
beta_info = beta_info(non_baseline_rows);

% exclude outliers

% z-scores each voxel (column)
beta_z = raj_fast_zscore(beta_mat(non_baseline_rows, :));
outliers = abs(beta_z) > z_score_cutoff;

% replace each outlier with the (signed) maximum allowable z score
beta_z(outliers) = z_score_cutoff * sign(beta_z(outliers));

fprintf('Excluding %d outliers (%.1f%%)\n', sum(outliers(:)), 100*mean(outliers(:)));

% count excluded betas by trial and by voxel
n_outliers_by_trial = sum(outliers, 2);
n_outliers_by_voxel = sum(outliers, 1);

n_voxels = sum(mask_mat(:));
trials_include_mask = n_outliers_by_trial/n_voxels < trial_exclude_outlier_proportion;
beta_z = beta_z(trials_include_mask, :);
beta_info = beta_info(trials_include_mask);
n_trials_excluded = sum(~trials_include_mask);
if n_trials_excluded
    fprintf('Excluding %d trials for having more than %.0f%% outliers\n', ...
        n_trials_excluded, round(trial_exclude_outlier_proportion*100));
end

%% cache these results in case we're called again with the same input
last_beta_z = beta_z;
last_beta_info = beta_info;
last_mask_mat = mask_mat;
last_outliers = outliers;
last_n_outliers_by_trial = n_outliers_by_trial;
last_n_outliers_by_voxel = n_outliers_by_voxel;
last_trials_include_mask = trials_include_mask;
last_subject = subject;
last_z_score_cutoff = z_score_cutoff;
last_trial_exclude_outlier_proportion = trial_exclude_outlier_proportion;