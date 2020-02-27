function vectorAnalysisPermTest(sessions, mode, Nperm, NstableVox, space, wordModel)
% vector analysis permutation test
% now only applied to whole-brain tests

% load permutation results
permResult_dir_template = fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
        'vectorAnalysis');
permResult_all = zeros(length(sessions), Nperm, 4);
for sess_idx = 1:length(sessions)
    session = sessions(sess_idx);
    permResult_dir = sprintf(permResult_dir_template, session.SubjectNum);
    load(fullfile(permResult_dir, ...
        sprintf('permResult_%s_%s_stable%d_perm%d_subj%d.mat', ...
        wordModel, mode, NstableVox, Nperm, sess_idx)), 'permResult');
    permResult_all(sess_idx, :, :) = permResult;
end

% random sample Nperm_all times from each fake subject to create the null
% distribution
Nperm_all = 10000;
permdist_all = zeros(4, Nperm_all);
for rank_idx = 1:4
    for i = 1:Nperm_all
        tmp_perm_rank = zeros(1, length(sessions));
        for subj = 1:length(sessions)
            % pick one fake subjects from the pool
            %tmp_perm_rank(subj) = permResult_all(subj, randperm(Nperm, 1), rank_idx);
            % use the fake subject one by one
            tmp_perm_rank(subj) = permResult_all(subj, i, rank_idx);
        end
        permdist_all(rank_idx, i) = mean(tmp_perm_rank);
    end
end
% load real ranking values
for sess_idx = 1:length(sessions)
    % load results
    session = sessions(sess_idx);
        load(sprintf(fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
            'vectorAnalysis', 'results_%s_%s_stable%d_perm%d.mat'), ...
            session.SubjectNum, wordModel, mode, NstableVox, Nperm));
   
    for i = 1:length(rankResult)
        for j = 1:length(rankResult{i})
            subjResult(j).subj_rank(i, sess_idx) = mean(mean(rankResult{i}(j).meanRelativeRank));
            subjResult(j).modelName= rankResult{i}(j).modelName;
        end
    end
end

% calculate mean
for i = 1:length(subjResult)
    subjResult(i).mean = mean(subjResult(i).subj_rank, 2);
    subjResult(i).std = std(subjResult(i).subj_rank, 0, 2);
    [h, subjResult(i).p] = ttest(subjResult(i).subj_rank, 0.5, 'Tail', 'right');
end

for rank_idx = 1:4
fprintf('permutation test %d p value is %f\n', rank_idx, ...
    sum(subjResult(rank_idx).mean < permdist_all(rank_idx, :))/Nperm_all);
end

% calculate within-subject p value
for rank_idx = 1:4
    numSigSubj = 0;
    for subj = 1:length(sessions)
        subjResult(rank_idx).withSubj_perm_p(subj) = sum(subjResult(rank_idx).subj_rank(subj) ...
            < permResult_all(subj, :, rank_idx)) / Nperm;
        if subjResult(rank_idx).withSubj_perm_p(subj) < 0.05
            numSigSubj = numSigSubj + 1;
        end
    end
    subjResult(rank_idx).numSigSubj = numSigSubj;
    
end
a=1;

