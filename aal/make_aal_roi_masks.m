function [roi_masks_struct] = make_aal_roi_masks(aal_dir, aal_filename, roi_subdir)
% function [roi_masks] = make_harvard_oxford_roi_masks(HarvOxf_dir, HarvOxf_filename, roi_subdir)
% reads in a harvard-oxford atlas image that has been rescaled into a particular
% voxel space, and creates new images that have cleaned up ROIs, separately for
% left and right hemispheres, and saves an image that has all the ROIs from both
% hemispheres.
%
% Input:
%   HarvOxf_dir: root directory of atlas
%   HarvOxf_filename: filename of specific atlas image to extract ROIs from
%   [roi_subdir]: name of subdirectory to save ROI images to.  If relative path,
%     will be relative to HarvOxf_dir.  Optional, default is no subdirectory.
%   isCort: 1 for cortical atlas, 0 for subcortical atlas
% Output:
%   [roi_masks_struct]: if asked for, will return a struct array, with fields
%   left, right, and both, containing logical masks for the corresponding ROI.
%
% Written by Rajeev Raizada (??/??/??)
% Modified by Dave Kleinschmidt
%   11/13/2014: Functionify

if nargin < 1
    aal_dir = fullfile(BaseDir(), 'aal');
    warning('Using default directory.  Are you sure this is what you want?\n\n%s', ...
        aal_dir);
end

starting_dir = pwd();
cleanup_object = onCleanup(@() cd(starting_dir));
cd(aal_dir);

spm('Defaults', 'fmri');


if nargin < 2
    error('Need to specify which particular aal atlas filename to use');
end

% defult to write ROI images in base HO atlas dir.
if nargin < 3
    roi_subdir = '';
end


if ~isempty(roi_subdir) && ~exist(roi_subdir, 'dir')
    fprintf('Creating subdirectory for ROI images');
    mkdir(roi_subdir);
end

V_HarvOxf = spm_vol(fullfile(aal_dir, aal_filename));
m_HarvOxf = spm_read_vols(V_HarvOxf);
[~, roi_prefix] = fileparts(aal_filename);

% little function to consistently create ROI filenames
make_roi_filename = @(leftrightall, roi_num) ...
    sprintf(fullfile(roi_subdir, [roi_prefix, '_%s_roi_%d.nii']), leftrightall, roi_num);

fprintf(['Will write images in files with names like\n\n  ' ...
         make_roi_filename('<left/right>', 1) '\n']);


roi_names = read_aal_names();
roi_masks_struct = struct();

for roi_num = 1:length(roi_names)

    fprintf('ROI num: % 2d (%s)\n', roi_num, roi_names(roi_num).name);
    this_roi = (m_HarvOxf==roi_names(roi_num).index);

    %%% Get rid of straggling disconnected voxels in the wrong place
    %[connected_region_labels,num_connected_regions] = spm_bwlabel(double(this_roi),18);
    %[freqs,bin_centers] = hist(nonzeros(connected_region_labels),[1:300]);
    %%% There are at most two proper regions, but sometimes only one if it
    %%% crosses the midline.
    %%% We will define a cutoff point of the number of voxels
    %voxel_count_cut_off = max(freqs)/2;
    %real_region_labels = find( freqs > voxel_count_cut_off );

    %this_roi_cleaned_up = ismember(connected_region_labels,real_region_labels);

    %leftV_this_roi = V_HarvOxf;
    %leftV_this_roi.fname = make_roi_filename('left', roi_num);
    %spm_write_vol(leftV_this_roi,left_roi);

    %rightV_this_roi = V_HarvOxf;
    %rightV_this_roi.fname = make_roi_filename('right', roi_num);
    %spm_write_vol(rightV_this_roi,right_roi);
    if strcmp(roi_names(roi_num).left_right, 'L')
        roi_masks_struct(ceil(roi_num/2)).left = logical(this_roi);
    else
        roi_masks_struct(ceil(roi_num/2)).right = logical(this_roi);
    end
    
    roi_masks_struct(ceil(roi_num/2)).name = roi_names(roi_num).name;
end;

for roi_num = 1:length(roi_masks_struct)
    roi_masks_struct(roi_num).both = roi_masks_struct(roi_num).left | roi_masks_struct(roi_num).right;
end
%%%% Make a volume containing all the ROIs
%{
V_all = V_HarvOxf;
V_all.fname = fullfile(roi_subdir, [roi_prefix, '_leftright_all.nii']);
m_all = zeros(size(m_HarvOxf));

for roi_num = 1:num_HarvOxf_rois,
    leftV_this_roi.fname = make_roi_filename('left', roi_num);
    leftV_this_roi = spm_vol(leftV_this_roi.fname);
    m_left = spm_read_vols(leftV_this_roi);

    rightV_this_roi.fname = make_roi_filename('right', roi_num);
    rightV_this_roi = spm_vol(rightV_this_roi.fname);
    m_right = spm_read_vols(rightV_this_roi);

    m_all = m_all + roi_num*m_left + (num_HarvOxf_rois+roi_num)*m_right;
end;
spm_write_vol(V_all,m_all);
%}
save(fullfile(roi_subdir, [roi_prefix, '_all_rois.mat']), 'roi_masks_struct');