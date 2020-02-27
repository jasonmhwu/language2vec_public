function batch_contrasts(subjects, analysis_subdir)

base_dir = BaseDir();
%%% Where does the SPM analysis happen? This is where teh SPM.mat file will
%%% be written, and is a subdirectory of the session's directory    
if nargin < 2
    analysis_subdir = 'SPM_analysis/';
end
% loop through subjects
for this_subject = subjects'
    
    % task 1: get to the right place.. directory wrangling tiimmmmee!!
    
    %%% Template to get a subject's analysis base directory:
    session_base_dir_template = fullfile(base_dir(), 'subjects', 'Subject%03d');
    % evaluate template by sprintf(session_base_dir_template, subnum, prepost)
    
    % put it all together and what do you get???
    analysis_dir_absolute = fullfile(...
        sprintf(session_base_dir_template, this_subject.SubjectNum), ...
        analysis_subdir);
    
    % and now.. contrasts! oohhh yyyeeaaahhhh!
    contrasts = make_contrasts(this_subject);
    config = this_subject;
    config.contrasts = contrasts;
    
    batch_obj = job_contrasts(config, analysis_dir_absolute);
    
    % save to a .mat file for record keeping purposes
    batch_obj_filename = fullfile(analysis_dir_absolute, 'batch_contrasts.mat');
    save(batch_obj_filename, 'batch_obj');
    
    % run the job
    try
        fprintf('\n\n================================================================================\n');
        fprintf('Running batch contrasts on subject %d\n', this_subject.SubjectNum);
        spm_jobman('initcfg');
        spm_jobman('serial', batch_obj);
    catch exception
        fprintf('Encountered a problem:\n%s\n  Skipping this subject', ...
            getReport(exception));
    end
    
end