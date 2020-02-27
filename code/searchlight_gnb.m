function [searchlight_test_class, searchlight_accuracies, searchlight_test_dist] = ...
         searchlight_gnb(activation_matrix, ...
                         training_times, testing_times, ...
                         searchlights, ...
                         ground_truth)
% function [searchlight_test_class, searchlight_accuracies, searchlight_test_dist] = ...
%          searchlight_gnb(activation_matrix, ... 
%                          training_times, testing_times, ...
%                          searchlight_indices, ...
%                          ground_truth)
% 
% Perform a Gaussian Naive Bayesian classifier analysis on activation data
% at specified training and testing times (rows) and broken out by
% specified searchlights (columns).
%
% Input: 
%   activation_matrix: a time/beta/condition x voxel activation matrix.
%     Rows are the different times, betas, or conditions to classify, and
%     columns are voxels
%   training_times: a cell vector which has the times to train each class
%     of the classifier on. That is, condition 1 training times should be
%     specified in training_times{1}, condition 2 in training_times{2}, etc.
%     Times are row indices in the activation_matrix
%   testing_times: a vector of the times (row indices) to be classified
%     into one of the classes from training_times.
%   searchlight_indices: a cell vector of searchlights.  Each searchlight
%     is a vector of the indices it contains (column indices in
%     activation_matrix).  This is the format that is output by
%     make_searchlights()
%   ground_truth: a vector of the actual classifications of the testing times.
%     Values should be the index of the correct training class (e.g. =1 for same
%     category as training_times{1}, =2 for trainign_times{2}, etc.)
%
% Output: 
%   searchlight_test_class: a testing_times x searchlight matrix of
%     test trial classifications.  Will be =1 for same as training_times{1}, =2
%     for training_times{2}, etc.
%   searchlight_accuracies: if a ground_truth is provided, then the overall
%     accuracy will be reported for each searchlight in a vector.

% Validate input: 

% is training times a length-2 cell? Only doing two classes for backwards
% compatibility...
if ~ iscell(training_times) || length(training_times) ~= 2
    error('Currently only supports two conditions');
end

% if ground truth is specified, is it the right length? 
if nargin > 4 && length(ground_truth) ~= length(testing_times)
    error('Mismatching length for ground_truth (%d) and testing_times (%d)', ...
        length(ground_truth), length(testing_times));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Train" the classifier: 

fprintf('  Training classifier...');
% Pull out training times for conditions 1 and 2
cond1_training_times = training_times{1};
cond2_training_times = training_times{2};
% For GNB, we only need to calculate the mean and variance for each
% cond once for each rep. Then we can re-use them for all spheres
cond1_training_means = mean(activation_matrix(cond1_training_times,:));
cond1_training_stds = std(activation_matrix(cond1_training_times,:));
cond2_training_means = mean(activation_matrix(cond2_training_times,:));
cond2_training_stds = std(activation_matrix(cond2_training_times,:));

%%% If we are using a frequency prior for the classes,
%%% then calculate the relative frequencies here
using_frequency_prior = 0;
if using_frequency_prior,
   num_training_times = length(cond1_training_times) + length(cond2_training_times);
   cond1_prior_p_value = length(cond1_training_times) / num_training_times;
   cond2_prior_p_value = length(cond2_training_times) / num_training_times;
end;

fprintf('done\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get voxel-by-voxel likelihoods, which can be pooled later by simply
% multiplying the likelihoods from each voxel in the searchlight.  This means we
% only need to calculate the likelihood of the test stimuli under each class for
% each voxel once.  (This is a nice advantage of the Gaussian Naive Bayes
% classifier).

fprintf('  Calculating test likelihoods for each voxel...');
% Now we can calculate the likelihoods for each testing point belonging either
% to cond1 or cond2, for all voxels at once
num_testing_times = length(testing_times);
num_mask_voxels = length(searchlights);

% We'll store the likelihoods for all testing time-points for all voxels
cond1_zscore_lhoods = zeros(num_testing_times,num_mask_voxels);
cond2_zscore_lhoods = zeros(num_testing_times,num_mask_voxels);

for testing_time_point_idx = 1:num_testing_times,
   testing_time_point = testing_times(testing_time_point_idx);
   testing_brainpattern = activation_matrix(testing_time_point,:);
   %%% Get the likelihood of this testing point for each voxel.  This
   %%% essentially sees how far away the testing pattern is from condition mean
   %%% mattern, scaled by the standard deviation 'pattern' (i.e. a z-score), and
   %%% then calcualte the likelihood of observing that z-score under a standard
   %%% normal distribution for all voxels
   cond1_zscore_lhoods(testing_time_point_idx,:) = ...
      normpdf(testing_brainpattern,cond1_training_means,cond1_training_stds);
   cond2_zscore_lhoods(testing_time_point_idx,:) = ...
      normpdf(testing_brainpattern,cond2_training_means,cond2_training_stds);

end;  %%% End of loop through testing-points

% convert to log-likelihood, which is numerically more robust (no small number
% multiplication problems) and combines by addition rather than multiplication.
cond1_log_lhoods = log(cond1_zscore_lhoods);
cond2_log_lhoods = log(cond2_zscore_lhoods);

fprintf('done\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate output for each searchlight by pooling single-voxel likelihoods
% computed above.
%

fprintf('  Processing searchlights');

% Main output will be a testing time x searchlight matrix of classifications: 
searchlight_test_class = zeros(num_testing_times, num_mask_voxels);
searchlight_test_dist = zeros(num_testing_times, num_mask_voxels);

% If ground truth is provided, we'll need to store the testing accuracy in a
% vector:
searchlight_accuracies = zeros(1, num_mask_voxels);

for voxel_num = 1:num_mask_voxels,

   sphere_inds_for_this_voxel = searchlights{voxel_num};

   cond1_sphere_log_lhoods = cond1_log_lhoods(:,sphere_inds_for_this_voxel);
   cond2_sphere_log_lhoods = cond2_log_lhoods(:,sphere_inds_for_this_voxel);
   % For all time points, sum the log p-vals for the voxels in the sphere
   % (calculating the overall likelihood under each condition's distribution)
   cond1_sphere_lhoods_sum = sum(cond1_sphere_log_lhoods,2);
   cond2_sphere_lhoods_sum = sum(cond2_sphere_log_lhoods,2);

   % If using a prior, adjust log-likelihood by that to get log-posterior.
   if using_frequency_prior,
      cond1_sphere_lhoods_sum = cond1_sphere_lhoods_sum + log(cond1_prior_p_value);
      cond2_sphere_lhoods_sum = cond2_sphere_lhoods_sum + log(cond2_prior_p_value);
   end;

   % If the cond2 log-likelihood sum is > cond1 sum, then assign output to
   % cond2, which is coded as +2. Cond1 is coded as +1
   testing_Category_outputs = ( cond2_sphere_lhoods_sum > cond1_sphere_lhoods_sum )+1;
   testing_Dist = (cond2_sphere_lhoods_sum - cond1_sphere_lhoods_sum);
   
   % Store classification of each test item by this searchlight in output mat
   searchlight_test_class(:, voxel_num) = testing_Category_outputs;
   searchlight_test_dist(:, voxel_num) = testing_Dist;
   
   % actual correct classification is provided through (optional) input
   % argument ground_truth
   if nargin > 4
       hasTruth = ~isnan(ground_truth);
       testing_Prop_correct = mean(testing_Category_outputs(hasTruth)==ground_truth(hasTruth));
       % Store accuracy in output matrix
       searchlight_accuracies(voxel_num) = testing_Prop_correct;
   end
   
   if rem(voxel_num, floor(num_mask_voxels / 10))==0
       fprintf('.');
   end
%    
%    if rem(voxel_num,floor(num_mask_voxels/3))==0,
%       disp(['Subj ' num2str(subj_num) ' rep ' num2str(rep) ...
%             ' test-set correct for voxel ' num2str(voxel_num) ...
%             ' out of ' num2str(num_mask_voxels) ...
%             ... %' at [' num2str([center_x center_y center_z]) ']' ...
%             ' with radius ' num2str(sphere_radius) ...
%             ' = ' num2str(testing_Prop_correct) ]);
%    end;

end;  %%% End of loop through searchlights

fprintf('done\n');