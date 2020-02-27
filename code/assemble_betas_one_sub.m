function [all_betas_mat, all_betas_info, mask_mat] = assemble_betas_one_sub(subj_num, save_output, reuse_if_available, mask_beta_subdir)
% function [all_betas_mat, all_betas_info, mask_mat] = assemble_betas_one_sub(subj_num, save_output, reuse_if_available)
% 
% Reads in individual beta image files output by SPM, and converts them
% into a beta index (trial) x mask voxel matrix of betas for later analysis
% by classifier etc.  Second return argument is a cell vector of the SPM
% header structs for each beta image (row in the resulting matrix).  If
% specified (save_output) output will be saved to a .mat file for later use
% (defaults to not save).

% define the output variables as persistent, in addition to the last-loaded
% subject/prepost.
persistent last_all_betas_mat last_all_betas_info last_mask_mat last_subj_num last_mask_beta_subdir

base_dir = BaseDir();
session_base_dir_template = fullfile(base_dir, 'subjects', 'Subject%03d');

if nargin < 1
    error('subject number or pre vs post not specified, but required');
end

% default to NOT save output to file
if nargin < 2
    save_output = 0;
end

if nargin < 3
    reuse_if_available = 0;
end
if nargin < 4
    mask_beta_subdir = 'SPM_analysis_motionReg';
end
%% use last loaded output if matches current request
import functools.*
persistents_empty = map(@isempty, ...
                        {last_subj_num, last_all_betas_mat, ...
                         last_all_betas_info, last_mask_mat, last_mask_beta_subdir});

if ~any([persistents_empty{:}]) && ...
        last_subj_num == subj_num && ...
        strcmp(last_mask_beta_subdir,  mask_beta_subdir)
    all_betas_mat = last_all_betas_mat;
    all_betas_info = last_all_betas_info;
    mask_mat = last_mask_mat;
    mask_beta_subdir = last_mask_beta_subdir;
    return
end

session_dir = sprintf(session_base_dir_template, subj_num);
output_filename = fullfile(session_dir, 'masked_betas_matrix.mat');

%% short circuit if existing save file is available
if reuse_if_available && exist(output_filename, 'file')==2
    % check whether any of the beta images are new
    beta_files = dir(fullfile(session_dir, mask_beta_subdir, 'beta_*.nii'));
    output_file = dir(output_filename);
    num_newer_betas = sum([beta_files.datenum] > output_file.datenum);
    if num_newer_betas
        fprintf('Cached results file exists, but there are %d newer beta image files. \n  Regenerating...\n  ', num_newer_betas);
    else
        fprintf('Reading cached results from %s...', fileparts(output_filename));
        load(output_filename);
        fprintf('done\n');

        last_subj_num = subj_num;
        last_mask_beta_subdir = mask_beta_subdir;
        last_all_betas_mat = all_betas_mat;
        last_all_betas_info = all_betas_info;
        last_mask_mat = mask_mat;

        return
    end
end

%% read in mask image as logical
maskV = spm_vol(fullfile(session_dir, mask_beta_subdir, 'mask.nii'));
mask_mat = spm_read_vols(maskV);

% get a list of x,y,z triplets for voxels in the mask (to be passed to
% spm_get_data below).
[mask_x, mask_y, mask_z] = ind2sub(size(mask_mat), find(mask_mat));
mask_xyz = [mask_x, mask_y, mask_z]';

%% read in beta images
% spm_select returns a character matrix of filenames
beta_filenames = spm_select('FPList', fullfile(session_dir, mask_beta_subdir), 'beta_.*\.nii');

% spm_vol can take a character matrix of filenames like spm_select returns, and
% returns a struct array for all the files, which can then be passed to the
% spm_read_vols function. (20% speed advantage over manually looping).
all_betas_info = spm_vol(beta_filenames);

% parse filenames to include only trial images
beta_info_table = struct2table(all_betas_info);
IsTrialBeta = contains(beta_info_table.descrip, 'bf');
all_betas_info = all_betas_info(IsTrialBeta);
num_betas = size(all_betas_info, 1);
fprintf('Reading %d beta images...', num_betas);

% use spm_get_data to just read the voxels in the mask from each image
all_betas_mat = spm_get_data(all_betas_info, mask_xyz);
fprintf('done\n');

%% save the output in a .mat file
if save_output
    fprintf('  saving output...');
    save(output_filename, ...
        'all_betas_mat', 'all_betas_info', 'mask_mat');
    fprintf('\n');
end

%% record the inputs and outputs in the persistent variables to speed up next time
last_subj_num = subj_num;
last_all_betas_mat = all_betas_mat;
last_all_betas_info = all_betas_info;
last_mask_mat = mask_mat;
last_mask_beta_subdir = mask_beta_subdir;
