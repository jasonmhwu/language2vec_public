% identify possible good ROI for each task
rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
                                     'harvard_oxford_animal_space.nii');

sigROIs = cell(1, 4);
groupResults = cell(5, 4);
voxSizeList = [100, 200, 400, 800, 1600];
for vox = 1:length(voxSizeList)
    load(sprintf('group_results_none_harvard-oxford_beta_stable%d_perm2.mat', voxSizeList(vox)));
    for mdl = 1:4
        groupResults{vox, mdl} = find(subjResult(mdl).p < 0.05);
    end
end

for i = 1:4
    sigROIs{i} = groupResults{1, i};
    for j = 2:3
        sigROIs{i} = intersect(sigROIs{i}, groupResults{j});
    end
end


% is any ROI better than whole brain?
BetterROI = zeros(5, 4, 96);
JustSigROI = zeros(5, 4, 96);
for vox = 1:length(voxSizeList)
    % load whole brain
    load(sprintf('group_results_none_wholebrain_beta_stable%d_perm2.mat', voxSizeList(vox)));
    wholeBrain_p = zeros(1, 4);
    for i = 1:4
        wholeBrain_p(i) = subjResult(i).p;
    end
    wholeBrain_p
    % load harvard oxford
    load(sprintf('group_results_none_harvard-oxford_beta_stable%d_perm2.mat', voxSizeList(vox)));
    for mdl = 1:4
        BetterROI(vox, mdl, :) = (subjResult(mdl).p < wholeBrain_p(mdl) & (subjResult(mdl).p < 0.05));
        JustSigROI(vox, mdl, :) = ((subjResult(mdl).p < 0.05));
    end
end

%% draw Raj's plot
sessions = subj_info();
sessions = sessions(1:end-1);
mode = 'whole-brain';
wordModel = 'none';
NstableVox = 1600;
Nperm = 2;
space = 'beta';
figure('Position', [0, 0, 2000, 1000]);
for sess_idx = 1:length(sessions)
    % load results
    session = sessions(sess_idx);
    if ~strcmp(mode, 'whole-brain')
        load(sprintf(fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
            'vectorAnalysis', 'results_%s_%s_stable%d_perm%d.mat'), ...
            session.SubjectNum, wordModel, mode, NstableVox, Nperm));
    else % whole-brain
        load(sprintf(fullfile(BaseDir(), 'subjects', 'Subject%03d', ...
            'vectorAnalysis', 'results_%s_wholebrain_%s_stable%d_perm%d.mat'), ...
            session.SubjectNum, wordModel, space, NstableVox, Nperm));
    end
    subplot(4, 4, sess_idx);
    y1 = rankResult{1}(1).meanRelativeRank(:);
    y2 = rankResult{1}(3).meanRelativeRank(:);
    hold on;
    for i = 1:length(y1)
        line([1, 2], [y1(i), y2(i)]);
    end
    hold off;
    axis([0.5, 2.5, 0, 1]);
    xticks(1:2);
    xticklabels({'true word', 'average of same-category words'});
    ylabel('ranking position');
    xlabel(sprintf('subject %d', sess_idx));
a = 1;

end