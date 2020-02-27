function [betas_mat_zscore, runInfo, mask_mat] = ...
            retrieve_actualTrials_one_sub(session)
    % retrieve the beta images of actual trials, clean up run info, and
    % standardize image by each run
    
    [all_betas_mat, all_betas_info, mask_mat] = ...
            assemble_betas_one_sub(session.SubjectNum, 1, 1, 'SPM_analysis_motionReg');
    % standardize images by each run
    runInfo = parseAndCombineOneSession(session.SubjectNum);
    runIdx = getColumn(runInfo, 'runNum');
    actualTrials = getColumn(runInfo, 'actualTrials');
    runIdx_uniq = unique(runIdx);
    for i = 1:length(runIdx_uniq)
        trials = find(runIdx == runIdx_uniq(i) & actualTrials);
        tmp = all_betas_mat(trials, :);
        tmp = (tmp - mean(tmp)) ./ std(tmp);
        betas_mat_zscore(length(trials) * (i-1) + (1:length(trials)), :) = tmp;
    end
    
    % clean up runInfo
    for i = 1:length(runInfo.varnames)
        runInfo.data{i} = runInfo.data{i}(actualTrials);
    end