function [searchlight_test_class, searchlight_accuracies] = ...
         searchlight_multi_gnb(activation_matrix, ...
                         training_times, testing_times, ...
                         searchlights, ...
                         ground_truth, ...
                         varargin)
% function [searchlight_test_class, searchlight_accuracies] = ...
%          searchlight_multi_gnb(activation_matrix, ... 
%                          training_times, testing_times, ...
%                          searchlight_indices, ...
%                          ground_truth)
% 
% Multiclass Gaussian Naive Bayes classifier analysis on activation data
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
%   varargin: all other arguments are passed to svmtrain
%
% Output: 
%   searchlight_test_class: a testing_times x searchlight matrix of
%     test trial classifications.  Will be =1 for same as training_times{1}, =2
%     for training_times{2}, etc.
%   searchlight_accuracies: if a ground_truth is provided, then the overall
%     accuracy will be reported for each searchlight in a vector.

n_searchlights = length(searchlights);

% initialize output
searchlight_test_class = zeros(length(testing_times), n_searchlights);
searchlight_accuracies = zeros(n_searchlights, 1);

% convert to column vectors
training_times = functools.map(@(x) x(:), training_times);

% pull out the training rows of the activation matrix
training_matrix = activation_matrix(cat(1, training_times{:}), :);
% create ground-truth labels
temp = cell(size(training_times));
for i = 1:length(training_times)
    temp{i} = i*ones(size(training_times{i}));
end
training_class = cat(1, temp{:});
% pull out testing rows of the activation matrix
testing_matrix = activation_matrix(testing_times, :);




% run classifier on every searchlight/roi:
parpool_exists = ~isempty(gcp('nocreate'));
if parpool_exists
    parfor sl_idx = 1:n_searchlights
            [this_class, this_acc] = ...
                do_one_searchlight(training_matrix(:, searchlights{sl_idx}), ...
                                   training_class, ...
                                   testing_matrix(:, searchlights{sl_idx}), ...
                                   ground_truth, ...
                                   varargin{:});
            searchlight_test_class(:, sl_idx) = this_class;
            searchlight_accuracies(sl_idx) = this_acc;
    end
else
    % set up logger
    str_temp = 'Searchlight % 6.d\n';
    str_clear = repmat('\b', 1, length(sprintf(str_temp, 0)));
    fprintf(str_temp, 0);
    for sl_idx = 1:n_searchlights
        % log every 10 searchlights
        if rem(sl_idx, 10) == 0
            fprintf([str_clear str_temp], sl_idx);
        end

        [searchlight_test_class(:, sl_idx), searchlight_accuracies(sl_idx)] = ...
            do_one_searchlight(training_matrix(:, searchlights{sl_idx}), ...
                               training_class, ...
                               testing_matrix(:, searchlights{sl_idx}), ...
                               ground_truth, ...
                               varargin{:});
    end
end


function [test_class, test_acc] = do_one_searchlight(sl_training_matrix, ...
                                                  training_class, ...
                                                  sl_testing_matrix, ...
                                                  ground_truth, ...
                                                  varargin)
% helper function to run on one searchlight

try 
    % Train classifier on this searchlight
    mgnb_mod = fitcnb(sl_training_matrix, training_class);
    test_class = predict(mgnb_mod, sl_testing_matrix);
    test_acc = mean(test_class == ground_truth);
catch
    test_class = nan(size(sl_testing_matrix,1), 1);
    test_acc = nan;
end