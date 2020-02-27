function stableVoxels(sessions)
% estimate stable Voxels for each wordType
% should do standardization first

if nargin < 1
    sessions = subj_info();
end
output_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'stableVoxels');
for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    
    [all_betas_mat, runInfo, mask_mat] = ...
            retrieve_actualTrials_one_sub(session);
    
    output_dir = sprintf(output_dir_template, session.SubjectNum);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    Stimulus = getColumn(runInfo, 'Stimulus');
    runIdx = getColumn(runInfo, 'runNum');
    Stim_uniq = unique(Stimulus);
    Nvox = size(all_betas_mat, 2);
    Nruns = max(runIdx);
    betaSumXRun = zeros(length(Stim_uniq), Nruns, Nvox);
    stability = zeros(1, Nvox);
    % metric: correlation sum across 6 runs
    for stim = 1:length(Stim_uniq)
        for run = 1:Nruns
            trial = runIdx == run & strcmp(Stimulus, Stim_uniq{stim});
            if sum(trial) ~= 1
                error('cannot find correct mapping');
            end
            betaSumXRun(stim, run, :) = all_betas_mat(trial, :);
        end
    end
    for vox = 1:Nvox
        corrXRun = triu(corr(betaSumXRun(:,:,vox)), 1);
        stability(vox) = sum(corrXRun(:));
    end
    output_fname = 'corrXRuns.nii';
    write_brain(stability, mask_mat, [], fullfile(output_dir, output_fname));
    % metric: sum of standard deviation for each word
    %{
    for class_idx = 1:length(wordType)
        stimInClass = Stimulus(wordType{class_idx});
        stimInClass_uniq = unique(stimInClass);
        assert(length(stimInClass_uniq) == 19); % sanity check that only works now!!!
        
        sumStdEachWord = zeros(1, size(all_betas_mat, 2));
        % sum up the std's for each unique word, and do things
        for w = 1:length(stimInClass_uniq)
            betas = all_betas_mat(find(strcmp(Stimulus, stimInClass_uniq{w})), :);
            sumStdEachWord = sumStdEachWord + std(betas);
        end
        output_fname = ['stableVoxels_', wordType_name{class_idx}, '.nii'];
        write_brain(sumStdEachWord, mask_mat, [], fullfile(output_dir, output_fname));
        output_neg_fname = ['stableVoxels_neg_', wordType_name{class_idx}, '.nii'];
        write_brain(-sumStdEachWord, mask_mat, [], fullfile(output_dir, output_neg_fname));
        
    end
    %}
end