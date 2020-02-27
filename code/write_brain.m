function [outV] = write_brain(voxels, mask, templateV, filename, descrip)
% write_brain(voxels, mask, templateV, filename[, descrip])
%
% writes the vector of voxel values, wrapped in mask, as a .nii file with
% provided filename and optionally the given description in the .descrip
% field.
%
% returns the SPM nifti header struct
% 
% to write a VECTOR of voxels, plugged into the logical-true values in mask,
%     write_brain(voxel_vector, mask, header, filename[, descrip])
% 
% to write a 3D volume, mask is ignored and can be empty:
%     write_brain(voxel_array, [], header, filename[, descrip])
    
if nargin < 5
    descrip = '';
end

if isempty(templateV)
    templateV = spm_vol(fullfile(BaseDir(), 'harvard_oxford',...
                                 'harvard_oxford_animal_space.nii'));
end

out_dir = fileparts(filename);
if ~ exist(out_dir, 'dir')
    mkdir(out_dir);
end

if prod(size(voxels)) == length(voxels)
    % vector of voxels:

    n_voxels = length(voxels);
    n_mask = sum(mask(:));

    % check that sizes are correct
    if n_voxels ~= n_mask
        error('Wrong number of voxels (%d) for mask (%d)', n_voxels, n_mask);
    end

    if ~ islogical(mask)
        mask = logical(mask);
    end

    outBrain = nan(size(mask));
    outBrain(mask) = voxels;
    
else 
    % three-d volume?
    if ~ all(size(voxels) == templateV.dim)
        error('Size mismatch between header and image');
    end
    
    outBrain = voxels;
end

outV = templateV;
outV.dt = [16 0];
outV.fname = filename;
outV.descrip = descrip;

outV = spm_write_vol(outV, outBrain);

end