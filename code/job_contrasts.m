function [matlabbatch] = job_contrasts(config, analysis_dir_absolute)
% construct an SPM batch job object for specifying contrasts.  Config is a
% struct that includes
%
% .subject: subject ID ('a'..'m')
% .contrasts: struct with contrast info.
%     .name: name of contrat
%     .vec:  contrast vector

spm_mat_location = fullfile(analysis_dir_absolute, 'SPM.mat');

consess = cell(size(config.contrasts));



%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_location};

for cond_i = 1:length(consess)
    this_contrast = config.contrasts(cond_i);
    matlabbatch{1}.spm.stats.con.consess{cond_i}.tcon.name = this_contrast.name;
    matlabbatch{1}.spm.stats.con.consess{cond_i}.tcon.convec = this_contrast.vec;
    matlabbatch{1}.spm.stats.con.consess{cond_i}.tcon.sessrep = 'none';
end

matlabbatch{1}.spm.stats.con.delete = 1;