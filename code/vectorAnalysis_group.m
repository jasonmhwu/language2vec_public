function vectorAnalysis_group(sessions, mode, Nperm, NstableVox, space, wordModel)
% can simple addition be a good model for specific areas?
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

if nargin < 1
    sessions = subj_info();
end
if nargin < 2  % 'harvard-oxford' or 'whole-brain'
    mode = 'harvard-oxford';
end
if nargin < 3
    Nperm = 1000;
end
if nargin <4
    NstableVox = 100;
end
if nargin < 5 % 'beta' or 'similarity' or 'simDiff'
    space = 'beta';
end
if nargin < 6
    wordModel = 'none';
end

[lr_rois, names] = load_rois(mode);

output_dir = fullfile(BaseDir(), 'group_analysis');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
subjResult(1).modelName = 'rankAllWords';
subjResult(2).modelName = 'rankWithinType';
subjResult(3).modelName = 'taxonimic rank';
subjResult(4).modelName = 'thematic rank';
for i = 1:4
    subjResult(i).subj_rank = zeros(length(lr_rois), length(sessions));
end

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
end

% do t-test for individual rois
for roi_idx = 1:length(lr_rois)
    for i = 1:length(subjResult)
        [~, subjResult(i).p(roi_idx)] = ttest(subjResult(i).subj_rank(roi_idx, :), ...
            0.5, 'Tail', 'right');
    end
end

% FDR correction across 96 rois
[h, crit_p, adj_ci_cvrg, subjResult(1).adj_p]=fdr_bh(subjResult(1).p);
[h, crit_p, adj_ci_cvrg, subjResult(2).adj_p]=fdr_bh(subjResult(2).p);
[h, crit_p, adj_ci_cvrg, subjResult(3).adj_p]=fdr_bh(subjResult(3).p);
[h, crit_p, adj_ci_cvrg, subjResult(4).adj_p]=fdr_bh(subjResult(4).p);


mdl_labels = {'Predict-identity', 'Predict-category', 'Predict-theme'};
mdl_colors = [255, 38, 0;
              95, 55, 255;
              0, 143, 0] / 255;
% figure 4: whole-brain 800 voxel analysis
if strcmp(mode, 'whole-brain')
    labelX = categorical(mdl_labels);
    labelX = reordercats(labelX, mdl_labels);
    y = [subjResult(1).mean, subjResult(3).mean, subjResult(4).mean];
    fig4 = bar(labelX, y);
    fig4.FaceColor = 'flat';
    fig4.CData(1, :) = mdl_colors(1, :);
    fig4.CData(2, :) = mdl_colors(2, :);
    fig4.CData(3, :) = mdl_colors(3, :);
    ylim([0.45, 0.6]);
    ylabel('Mean Relative Rank', 'FontSize', 15);
    xlabel('Ranking Metrics', 'FontSize', 15);
    title('Word predicted from addition and subtraction of fMRI patterns');
end

% figure 5: 8-ROI 100 voxel analysis
if strcmp(mode, 'general-semantic-network')
    labels = {'IPL', 'ITG&MTG', 'PHG', 'dmPFC', 'IFG', 'vmPFC', 'PCG', 'ATL'};
    labels_rearranged = {'IPL', 'ITG&MTG', 'PHG', 'PCG', 'ATL', 'IFG', 'dmPFC', 'vmPFC'};
    labelX = categorical(labels);
    labelX = reordercats(labelX, labels_rearranged);
    mdlList = [1, 3, 4];
    for mdl = 1:3
        subplot(3, 1, mdl);
        y = subjResult(mdlList(mdl)).mean;
        bar(labelX, y, 'FaceColor', [0.5, 0.5, 0.5]);
        ylim([0.45, 0.6]);
        yline(0.5, '-.k', 'LineWidth', 2);
        ylabel('Rank', 'FontSize', 15);
        title(mdl_labels(mdl), 'Color', mdl_colors(mdl, :));
    end
    xlabel('Semantic ROIs', 'FontSize', 15);
end

% write brain maps
if ~strcmp(mode, 'whole-brain')
    for mdl_idx = 1:length(subjResult)
        p_mask = fill_mask_voxels_with_roi_values(lr_rois, -log(subjResult(mdl_idx).p));
        write_brain(p_mask, [], [], fullfile(output_dir, ...
            sprintf('pVal_%s_stable%d_perm%d.nii', subjResult(mdl_idx).modelName, NstableVox, Nperm)));
    end
    save(sprintf(fullfile(BaseDir(), 'group_analysis', ...
        'group_results_%s_%s_%s_stable%d_perm%d.mat'), ...
        wordModel, mode, space, NstableVox, Nperm), 'subjResult');
else
    save(sprintf(fullfile(BaseDir(), 'group_analysis', ...
        'group_results_%s_whole-brain_%s_stable%d_perm%d.mat'), ...
        wordModel, space, NstableVox, Nperm), 'subjResult');
end
