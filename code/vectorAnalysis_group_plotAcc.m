function vectorAnalysis_group_plotAcc(mode, Nperm, voxSizeValues, space)
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

includeWord2Vec = 0;

lr_rois = load_rois(mode);

mdlSize = 4;
Nrois = length(lr_rois);
modelLegendName = {'true word vs. all', 'true word vs. same class', 'taxonomic rank', 'thematic rank'};
if includeWord2Vec
    Nrepre = length(voxSizeValues)*2 + 1;
    xLabelName = [string(voxSizeValues), "w2v", string(voxSizeValues)+"+w2v"];
else
    Nrepre = length(voxSizeValues);
    xLabelName = string(voxSizeValues);
end
p_value = zeros(mdlSize, Nrepre, Nrois);
means = zeros(mdlSize,  Nrepre, Nrois);
stds = zeros(mdlSize, Nrepre, Nrois);
subjRank = cell(mdlSize, Nrepre, Nrois);
for i = 1:length(voxSizeValues)
    load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'group_results_none_%s_%s_stable%d_perm%d.mat'), ...
            mode, space, voxSizeValues(i), Nperm), 'subjResult');
        for j = 1:mdlSize
            for roi = 1:Nrois
                [~, p_value(j, i, roi)] = ttest(subjResult(j).subj_rank(roi, :), 0.5, 'Tail', 'right');
                means(j, i, roi) = mean(subjResult(j).subj_rank(roi, :), 2);
                stds(j, i) = std(subjResult(j).subj_rank(roi, :), 0, 2);
                subjRank{j, i, roi} = subjResult(j).subj_rank(roi, :);
            end
        end
    if includeWord2Vec
        load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'group_results_word2vec_%s_%s_stable%d_perm%d.mat'), ...
            mode, space, voxSizeValues(i), Nperm), 'subjResult');
       for j = 1:mdlSize
            for roi = 1:Nrois
                [~, p_value(j, i+length(voxSizeValues)+1, roi)] = ttest(subjResult(j).subj_rank(roi, :), 0.5, 'Tail', 'right');
                means(j, i+length(voxSizeValues)+1, roi) = mean(subjResult(j).subj_rank(roi, :), 2);
                stds(j, i+length(voxSizeValues)+1) = std(subjResult(j).subj_rank(roi, :), 0, 2);
                subjRank{j, i+length(voxSizeValues)+1, roi} = subjResult(j).subj_rank(roi, :);
            end
       end
    end
end

if includeWord2Vec
load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'group_results_word2vec_%s_%s_stable0_perm%d.mat'), ...
            mode, space, Nperm), 'subjResult');
        for j = 1:mdlSize
            for roi = 1:Nrois
                [~, p_value(j, length(voxSizeValues)+1, roi)] = ttest(subjResult(j).subj_rank(roi, :), 0.5, 'Tail', 'right');
                means(j, length(voxSizeValues)+1, roi) = mean(subjResult(j).subj_rank(roi, :), 2);
                stds(j, length(voxSizeValues)+1) = std(subjResult(j).subj_rank(roi, :), 0, 2);
                subjRank{j, length(voxSizeValues)+1, roi} = subjResult(j).subj_rank(roi, :);
            end
        end
end
sems = stds ./ sqrt(length(subjResult(1).subj_rank));
% calculate t-test
ts = (means - 0.5) ./ (stds / sqrt(13));

% figure 1: mean relative rank vs. # voxels, all models
figure('Position', [0, 0, 1000, 800]);
for i = 1:mdlSize
    errorbar((1:size(means, 2)) + i * 0.07, means(i, :), sems(i, :));
    hold on;
end
title('Mean Relative Rank vs. number of stable voxels')
xlabel('number of stable voxels');
ylabel('mean relative rank');
xticks(1:length(xLabelName));
xticklabels(xLabelName);
legend(modelLegendName);
saveas(gca, ...
             fullfile(BaseDir(), 'figures', ...
                      'RankVoxelCurve_all_model.png'), 'png');
         
% figure 2: p value vs. # voxels, all models
figure('Position', [0, 0, 1000, 800]);
plot(p_value');
title('p value vs. number of stable voxels')
xlabel('number of stable voxels');
xticks(1:length(xLabelName));
xticklabels(xLabelName);
legend({'true word vs. all', 'true word vs. same class', 'taxonomic rank', 'thematic rank'});
saveas(gca, ...
             fullfile(BaseDir(), 'figures', ...
                      'PValueVoxelCurve_all_model.png'), 'png');

% figure 3-5: mean relative rank vs. # voxel, each model
% with individual subject point to see grouping
for mdl = 1:mdlSize
    figure('Position', [0, 0, 1000, 800]);
    plot(means(mdl, :));
    hold on;
    for i = 1:size(means, 2)
        plot(i * ones(1, length(subjRank{mdl, i})), subjRank{mdl, i}, 'o');
    end
    hold off;
    xlabel('number of stable voxels');
    ylabel('mean relative rank');
    xticks(1:length(xLabelName));
    xticklabels(xLabelName);
    title(sprintf('%s model: rank vs. number of voxels', modelLegendName{mdl}));
    saveas(gca, ...
             fullfile(BaseDir(), 'figures', ...
                      sprintf('RankVoxelCurve_%s.png', ...
                      modelLegendName{mdl})), 'png');

end


% figure 4: correlation plot between catch trial acc and model 1
% performance
subjects = subj_info(1);
subjects = subjects(1:end-1);
catchTrialAcc = zeros(1, length(subjects));
for subj = 1:length(subjects)
    subject = subjects(subj);
    catchTrialAcc(subj) = double(subject.Correct) / 18;
end
figure();
scatter(catchTrialAcc, subjRank{1, 4});
model = fitlm(catchTrialAcc, subjRank{1, 4});
hold on;
plot(catchTrialAcc, model.Coefficients.Estimate(1) + ...
    model.Coefficients.Estimate(2) * catchTrialAcc);
xlabel('Catch Trial Accuracy');
ylabel('True Word Against All Model Rank (800 voxels)');
title('Quality Check: Catch Trial Acc Does Not Predcit Model Rank');
text(0.9, 0.56, sprintf('slope: %f', model.Coefficients.Estimate(2)));
saveas(gca, ...
             fullfile(BaseDir(), 'figures', ...
                      sprintf('CatchAccRankCurve_%s.png', ...
                      modelLegendName{1})), 'png');

% figure 5: paired t-test between model1 & 3, mode 1 & 4
% model 1 is marignally significant against model 3
% model 4 is significant against model 4
pairedT = zeros(2, size(subjRank, 2));
for i = 1:size(pairedT, 2)
    [h, pairedT(1, i)] = ttest(subjRank{1, i}, subjRank{3, i}, 'Tail', 'right');
    [h, pairedT(2, i)] = ttest(subjRank{1, i}, subjRank{4, i}, 'Tail', 'right');
end

%%figure 6: paired plot between model 1 and 3 and 4
figure();
hold on;
for subj =1:13
    line([1, 2], [subjRank{1, 5}(subj), subjRank{3, 5}(subj)]);
    line([3, 4], [subjRank{1, 5}(subj), subjRank{4, 5}(subj)])
end
plot([0, 5], [0.5, 0.5], 'b:');
set(gca, 'FontSize', 20);
text(1.5, 0.565, sprintf('p = %.3f', pairedT(1, 5)), 'FontSize', 20);
text(3.5, 0.565, sprintf('* p = %.3f', pairedT(2, 5)), 'FontSize', 20);

xticks(1:4);
    xticklabels({'true word', 'Does predicted word have correct category', 'true word', 'Does predicted word have correct context'});
    ylabel('ranking position');
    title('for 1600 voxels in the whole brain, true word ranking is marginally higher than correct category ranking, and is significantly higher than correct context ranking');
