function vectorAnalysis_group_plotAcc_ROI(mode, Nperm, voxSizeValues, space)
% can simple addition be a good model for specific areas?
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

if nargin < 1  % 'harvard-oxford' or 'whole-brain'
    mode = 'harvard-oxford';
end
if nargin < 2
    Nperm = 1000;
end
if nargin <3
   voxSizeValues = 100;
end
if nargin < 4 % 'beta' or 'similarity' or 'simDiff'
    space = 'beta';
end

[lr_rois, names] = load_rois(mode);

mdlSize = 4;
Nrois = length(lr_rois);
Nsubjects = 13;
Nrepresentation = length(voxSizeValues);
modelLegendName = {'true word vs. all', 'true word vs. same class', 'taxonomic rank', 'thematic rank'};
xLabelName = [string(voxSizeValues)];
p_value = zeros(mdlSize, length(voxSizeValues), Nrois);
adj_p = zeros(mdlSize, length(voxSizeValues), Nrois);
means = zeros(mdlSize, length(voxSizeValues), Nrois);
ts = zeros(mdlSize, length(voxSizeValues), Nrois);
stds = zeros(mdlSize, length(voxSizeValues), Nrois);
subjRank = cell(mdlSize, length(voxSizeValues), Nrois);
for i = 1:length(voxSizeValues)
    load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'group_results_none_%s_%s_stable%d_perm%d.mat'), ...
            mode, space, voxSizeValues(i), Nperm), 'subjResult');
        for j = 1:mdlSize
            for roi = 1:Nrois
                [~, p_value(j, i, roi)] = ttest(subjResult(j).subj_rank(roi, :), 0.5, 'Tail', 'right');
                means(j, i, roi) = mean(subjResult(j).subj_rank(roi, :), 2);
                stds(j, i, roi) = std(subjResult(j).subj_rank(roi, :), 0, 2);
                subjRank{j, i, roi} = subjResult(j).subj_rank(roi, :);
            end
             [h, crit_p, adj_ci_cvrg, adj_p(j, i, :)]=fdr_bh(p_value(j, i, :), 0.05);
        end
end

repre = size(means, 2);
for repre = 1:Nrepresentation
    figure('Position', [0, 0, 1000, 1000]);
    ax = gca;
    hold on;
    
    % do significance test
    y = squeeze(means(:, repre, :))';
    y_std = squeeze(stds(:, repre, :))';
    
    % start plotting error bars
    ngroups = size(y, 1);
    nbars = size(y, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    p1 = bar(y, 'BarWidth', 1);
    
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
        errorbar(x, y(:,i), y_std(:, i)./sqrt(Nsubjects-1), 'k', 'linestyle', 'none');
    end
    % plot significant stars
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
        for j = 1:ngroups
            mysigstar(ax, x(j), ax.YLim(2)-0.05, p_value(i, repre, j));
        end
    end
    xticks(1:Nrois);
    xticklabels(names);
    ylabel("mean relative rank");
    legend(p1, modelLegendName);
    title(sprintf('Relative Rank for each ROI and each model, at %d voxels\n', voxSizeValues(repre)));
    saveas(gca, ...
             sprintf(fullfile(BaseDir(), 'figures', ...
                      'RankVoxelBar_%s_%d.png'), ...
                      mode, voxSizeValues(repre)), 'png');

end
a=1;

