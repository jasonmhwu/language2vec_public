function [sessions] = do_classifier(sessions, classifier_fcn, classes_fcn, job_name, ...
                                    searchlight_fcn, brain_filling_fcn, NstableVox, ...
                                    varargin)
% function [sessions] = do_classifier(sessions, classifier_fcn, classes_fcn, job_name
%                                     searchlight_fcn, brain_filling_fcn)
%
% Run the classifier analysis pipeline on the specified sessions.  Manages
% proper wrapping of input and saving output, and returns a sessions struct
% with information on where the output files are written.  
% 
% If you want to read in the full classifier structs, use
% 
% >> sessions = load_classifiers(sessions);
% 
% To read in the accuracy images as a four-D volume:
% 
% >> classifiers = [sessions.classifier];
% >> acc_vols_four_d = spm_read_vols([classifiers.acc_image_header]);
% 
% (or, if you have github.com/kleinschmidt/functools_for_matlab):
% 
% >> functools.pluck(sessions, 'classifier', 'acc_image_header')
% 
% Input: 
%   sessions: struct info of sessions (needs to have SubjectNum and PrePost)
%   classifier_fcn: function handle for classifier that is passed to
%     searchlight_classifier_crossval 
%   classes_fcn(session): function that takes one session and returns a cell
%     array with the trial/beta indices for each class.  See do_gnb_inout for
%     an example.
%   job_name: string with the job name.  will be added to the filenames of
%     any output and should uniquely identify the type of classifier analysis.
%   searchlight_fcn: how to generate subsets of voxels to classify.  Defaults to
%     generating searchlights (spherical region around each voxel), but any
%     function that takes a mask and returns a cell array of mask index vectors
%     will do.
%   brain_filling_fcn(rois, values, mask): how to go from a vector of values for
%     each searchlight to a brain image.  should return a vector of values
%     suitable for passing to write_brain(values, mask, ...);
%   beta_z_cutoff: if not missing or empty, beta values will be z-scored within
%     voxels, and any trial x voxel combination with absolute z-score greater
%     than the cutoff will be set to the (signed) cutoff value. For z-scoring
%     but no capping use Inf.
%   varargin: Additional arguments passed to classifier_fcn
% Output:
%   sessions: struct with one entry per input session, with field
%     "classifier" added with the filename of the saved classifier struct in
%     field .mat and the SPM Nifty volume header for the saved accuracy image
%     in .nii_header.
% 
% Side effects:
%   For each session, a .nii file is written with the searchlight accuracy
%   maps in LinNNN/<pre/post>/<job name>_LinNNN_<pre/post>.nii, and a .mat
%   file with the classifier struct is written in the same location with an
%   analogous filename.  The returned sessions struct is written to
%   cache/classifier/<job_name>.mat, and contains pointers to the accuracy
%   images in the classifier.acc_image_header fields.
    
import functools.*

fprintf('Running classifier analysis (job name: %s) on %d sessions\n\n', ...
        job_name, length(sessions));

if nargin < 5
    % default to 3-voxel radius searchlights
    sl_radius = 3;
    searchlight_fcn = partial(make_searchlights, sl_radius);
end

if nargin < 6
    brain_filling_fcn = [];
end

if nargin < 7 
    NstableVox = 50;
end

sessions = cell2mat(map(@(session) do_one(session, ...
                                          classifier_fcn, ...
                                          classes_fcn, ...
                                          job_name, ...
                                          searchlight_fcn, ...
                                          brain_filling_fcn, ...
                                          NstableVox, ...
                                          varargin{:}), ...
                        sessions));

fprintf('Saving final output sessions to cache/classifier/%s.mat...', ...
        job_name);

output_dir = fullfile(BaseDir(), 'cache', 'classifier');
if ~ exist(output_dir, 'dir')
    mkdir(output_dir);
end

save(fullfile(output_dir, sprintf('%s.mat', job_name)), 'sessions');
fprintf('done!\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%
% Run one session.
%%%%%%%%%%%%%%%%%%%%%%
function [session] = do_one(session, classifier_fcn, classes_fcn, job_name, ...
                            searchlight_fcn, brain_filling_fcn, ...
                            NstableVox, varargin) 

    import functools.*
    stableVoxel_fname_template = fullfile(BaseDir, 'subjects', 'Subject%03d', 'stableVoxels', 'corrXRuns.nii');

    fprintf('Session: Lin%03d\n\n', session.SubjectNum);
    
    % load in beta images
    [all_betas_mat, runInfo, mask_mat] = ...
        retrieve_actualTrials_one_sub(session);
    
    % temporary: shuffle the indices, to see if there are bugs
    %all_betas_mat = all_betas_mat(randperm(size(all_betas_mat, 1)), :);
   
    classifier.job_name = job_name;
    classifier.session_dir = make_session_dir(session);

    % use the classes_fcn to extract the true classes
    classifier.classes = classes_fcn(session);
    
    % concatenate together the list of trials in each class
    tmp = classifier.classes';
    classifier.trials_to_classify = cat(1, tmp{:});
    
    % create ground truth vector (1 for every trial in class 1, 2 for every
    % trial in class 2, etc.)
    classifier.ground_truth = [];
    for label = 1:size(classifier.classes, 1)
        for iden = 1:size(classifier.classes, 2)
            classifier.ground_truth = [classifier.ground_truth; ...
                label*ones(length(classifier.classes{label, iden}), 1)];
        end
    end
    
    % using leave-one-word-out cross-validation procedure
    classifier.n_cv_folds = size(classifier.classes, 2);
    
    % find the most stable voxels
    stable_fname = sprintf(stableVoxel_fname_template, session.SubjectNum);
    stableVoxels = spm_read_vols(spm_vol(stable_fname));
        
    % convert to index
    mask_inds = find(mask_mat);     
    num_mask_voxels = length(mask_inds);
    image_to_mask_indices = zeros(size(mask_mat));
    image_to_mask_indices(logical(mask_mat)) = 1:num_mask_voxels;
    % check that this works: 
    if ~ all(image_to_mask_indices(mask_inds)' == 1:num_mask_voxels)
        error('Image space to mask space inverse mapping is wrong somehow!');
    end
    
    if strcmp(searchlight_fcn, 'whole-brain')
        stabilityInRoi =sort(stableVoxels(logical(mask_mat)), 'descend');
        stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
        stableVox_mat = mask_mat & (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
       
        searchlights{1} = image_to_mask_indices(stableVox_mat);
        
    elseif strcmp(searchlight_fcn, 'harvard-oxford')
        load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
        lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
        
        stable_lr_rois = cell(1, length(lr_rois));
        for i = 1:length(lr_rois)
            stabilityInRoi =sort(stableVoxels(lr_rois{i}), 'descend');
            stabilityInRoi = stabilityInRoi(~isnan(stabilityInRoi));
            stable_lr_rois{i} = lr_rois{i} & ...
                (stableVoxels >= stabilityInRoi(min(NstableVox, length(stabilityInRoi))));
        end
        
        for i = 1:length(lr_rois)
            searchlights{i} = image_to_mask_indices( ...
                stable_lr_rois{i});
        end
    else % rad-3 searchlight
        searchlights = make_searchlights(3, mask_mat);
    end
    
    % run it.
    fprintf('Running classifier...\n');
    [classifier.cv_predict, classifier.cv_acc, classifier.partitions] = ...
        searchlight_classifier_crossval(all_betas_mat, ...
                                        classifier.trials_to_classify, ...
                                        classifier.n_cv_folds, ...
                                        searchlights, ...
                                        classifier.ground_truth, ...
                                        classifier_fcn, ...
                                        classifier.classes, ...
                                        varargin{:});
    
    fprintf('  Failed on %d%% of folds\n', ...
            round(100 * mean(isnan(classifier.cv_acc(:)))));
    
    % save output
    fprintf('Saving output for this session...');
    output_filename_stem = sprintf(fullfile(classifier.session_dir, ...
                                            '%s_Subject%03d_stable%d'), ...
                                   job_name, ...
                                   session.SubjectNum, ...
                                   NstableVox);
    
    % write searchlight accuracy maps (centered about 0.5) to a .nii file.
    % use one of the beta image headers as a template header...
    %template_header = all_betas_info(1);
    % description string to write in the .descrip field of the nifti header
    description = sprintf('classifier %s: Subject%03d', ...
                          job_name, session.SubjectNum);

    acc_img_voxels = mean(classifier.cv_acc) - 1/3;
    if ~isempty(brain_filling_fcn)
        acc_img_voxels = brain_filling_fcn(searchlights, ...
                                           acc_img_voxels, ...
                                           mask_mat);
    end
    
    try
        classifier.acc_image_header = ...
            write_brain(acc_img_voxels, mask_mat, [], ...
                        [output_filename_stem '.nii'], description);
        session.classifier.acc_image_header = classifier.acc_image_header;
    catch
        warning('Error while trying to write brain image. Omitting');
    end
    
    save([output_filename_stem '.mat'], '-struct', 'classifier');
    session.classifier.mat = [output_filename_stem '.mat'];
    
    fprintf('done!\n\n');
    
    % return session with classifier and acc image filenames/headers
    session.classifier.job_name = job_name;
end

function [session_dir] = make_session_dir(session)
    session_dir_template = fullfile(BaseDir(), ...
                                    'subjects', 'Subject%03d');
    session_dir = sprintf(session_dir_template, ...
                          session.SubjectNum);
end