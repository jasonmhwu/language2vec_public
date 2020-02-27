function [sphere_XYZ_indices_cell, mask_matrix] = make_searchlights(sphere_radius_range, mask_matrix, min_num_voxels)
% sphere_XYZ_indices_cell = make_searchlights(sphere_radius_range, mask_matrix, min_num_voxels)
% 
% Generate a cell array of searchlights for different radii.  Cell array is a
% mask_voxel x radius array, each entry of which is a vector with the indices of
% the MASK VOXELS which are included in that searchlight.  That is, the vector
% [1, 2, ...] refers to the first and second MASK VOXELS, rather than the first
% and second voxels in the whole image.  To convert to meatspace voxels, do
% somthing like
%
%   mask_inds(sphere_XYZ_indices_cell{1});
%
% The strategy used is to first calculate the relative indices (in image space)
% for a spotlight of each size requested.  This is done with sphere_indices()
% which calculates a sphere mask and then finds its indices in the space of the
% overall mask matrix, finally subtracting the index of the center voxel.  Then
% for each voxel in the mask, these relative indices are adjusted to find the
% searchlight sphere around that voxel in image space, which are then
% transformed into the indices of the voxels in the mask.
%
% The third and final argument is the minimum number of voxels that are allowed
% in each searchlight.  Any searchlights with less than the minimum are not
% returned.

% default to 3 voxel searchlight
if nargin < 1
    sphere_radius_range = 3;
end

% load the mask image if none is provided
if nargin < 2
    Vmask = spm_vol(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space.nii'));
    mask_matrix = spm_read_vols(Vmask);
end
% convert mask matrix to logical
mask_matrix = mask_matrix | 0; 

if nargin < 3
    min_num_voxels = 1;
end

%%% The indices in image space of where mask=1
mask_inds = find(mask_matrix);     
num_mask_voxels = length(mask_inds);

%%%% Create an array that has the image-to-mask (inverse) index map, so we
%%%% can convert from the overall image indices that are calculated by the
%%%% sphere_indices() function to the mask indices that this function
%%%% returns.
image_to_mask_indices = zeros(size(mask_matrix));
image_to_mask_indices(mask_matrix) = 1:num_mask_voxels;
% check that this works: 
if ~ all(image_to_mask_indices(mask_inds)' == 1:num_mask_voxels)
    error('Image space to mask space inverse mapping is wrong somehow!');
end

x_size = size(mask_matrix, 1); % Vmask.dim(1);
y_size = size(mask_matrix, 2); % Vmask.dim(2);
z_size = size(mask_matrix, 3); % Vmask.dim(3);

max_mask_index = x_size * y_size * z_size;

%%% pre-allocate results cell array, which is mask_voxel x radius
sphere_XYZ_indices_cell = cell(num_mask_voxels,length(sphere_radius_range));

%%% Pre-generate relative indices for all specified sphere radii
sphere_indices_cell = cellfun(@(r) make_sphere_indices(r, [x_size, y_size, z_size]), ...
                              num2cell(sphere_radius_range), ...
                              'UniformOutput', 0);


for voxel_num =1:length(mask_inds),

    if rem(voxel_num,10000)==0,
        disp([ 'Voxel number ' num2str(voxel_num) ' out of ' num2str(num_mask_voxels) ]);
    end;
    
    center_ind = mask_inds(voxel_num);
    
    for sphere_radius_idx = 1:length(sphere_radius_range),
        % adjust relative indices to center on current voxel
        sphere_inds = sphere_indices_cell{sphere_radius_idx} + center_ind;
        % remove voxels outside the mask array altogether
        sphere_inds = sphere_inds(sphere_inds > 0 & sphere_inds < max_mask_index);
        % pick out only the voxels that are included in the mask,
        % converting to mask voxel index (from image voxel index)
        sphere_XYZ_indices_cell{voxel_num, sphere_radius_idx} = ...
            image_to_mask_indices(sphere_inds(mask_matrix(sphere_inds)));
    end;
    
end;

% filter out searchlights with less than the minimum number of voxels
searchlight_sizes = cellfun(@length, sphere_XYZ_indices_cell);
if any(searchlight_sizes < min_num_voxels)
    warning('Removing %d searchlights with less than %d voxels', ...
            sum(searchlight_sizes < min_num_voxels), ...
            min_num_voxels);
    sphere_XYZ_indices_cell = sphere_XYZ_indices_cell(searchlight_sizes >= min_num_voxels);
    mask_matrix(mask_inds(searchlight_sizes < min_num_voxels)) = 0;
end

end

function [brain_space_sphere_inds_relative, sphere_x, sphere_y, sphere_z] = make_sphere_indices(radius, brain_space_size)
% sphere_inds = sphere_indices(radius, size)
%
% generate sphere of specified radius (in voxels), and return linear
% indices (and, optionally, xyz coordinates).  generates relative to a
% matrix of specified overall size (which determines the mapping from xyz to
% linear indices), and computes RELATIVE indices (where the index of the
% origin is 0)

%%

% radius rounded down because half a voxel isn't included in the sphere.
radius_int = floor(radius);
sphere_diameter = radius_int * 2 + 1;
center_coord = radius_int + 1;

% throw an error if the target matrix is too small to fit the whole sphere
% (we could just truncate but there's no clear way to do this best).
if any(brain_space_size < sphere_diameter)
    error('Target matrix size is smaller than sphere of radius');
end

%%% Raj verbose and slightly grating comments
%%% The code below makes a searchlight-sized block of distances.
%%% We will then slide
%% make 3d sphere mask and get indices of voxels in sphere
% start by generating a cube of the right size
coords = 1:sphere_diameter;
% ndgrid returns a 3d array that has the x-coordinates, one with
% y-coordinates, one with the z-coordinates
[x, y, z] = ndgrid(coords);
% then we can calculate the distance from the center for each voxel in
% the 3d array using element-wise operations
dist = sqrt((x-center_coord).^2 + (y-center_coord).^2 + (z-center_coord).^2);

% get the absolute indices of all the voxels in the cube within the
% specified radius
sphere_mask = dist <= radius;
sphere_inds_absolute = find(sphere_mask);

%% convert cube-based indices to destination indices, and get offset

% sphere_x = x(sphere_mask);
% sphere_y = y(sphere_mask);
% sphere_z = z(sphere_mask);

%%% An equivalent way of doing the above, using ind2sub
[sphere_x, sphere_y, sphere_z] = ind2sub(size(sphere_mask),sphere_inds_absolute);

% convert xyz to linear indices in destination matrix
% (this corresponds to the sphere nestled into the furthest corner of the
% brain space)
brain_space_sphere_inds_absolute = ...
    sub2ind(brain_space_size, sphere_x, sphere_y, sphere_z);

% convert to relative indices (offset from index of origin) unless
% specified otherwise
brain_space_sphere_inds_relative = ...
    brain_space_sphere_inds_absolute - sub2ind(brain_space_size, center_coord, center_coord, center_coord);

end