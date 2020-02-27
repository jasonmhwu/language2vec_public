function RDMAnalysis_group_plotAcc(mode, voxSizeValues, space)
% can simple addition be a good model for specific areas?
% first run stableVoxels analysis, then pick the most stable voxels for
% each ROI and perform vector Analysis

Nsubjects = 13;

if nargin < 1  % 'harvard-oxford' or 'whole-brain'
    mode = 'harvard-oxford';
end
if nargin <2
   voxSizeValues = 100;
end
if nargin < 3 % 'beta' or 'similarity' or 'simDiff'
    space = 'beta';
end

if strcmp(mode, 'harvard-oxford')
    %rois = make_harvard_oxford_roi_masks(fullfile(BaseDir(), 'harvard_oxford'), ...
    %                                 'harvard_oxford_animal_space.nii');
    load(fullfile(BaseDir(), 'harvard_oxford', 'harvard_oxford_animal_space_all_rois.mat'));
    lr_rois = {roi_masks_struct.left, roi_masks_struct.right};
elseif strcmp(mode, 'fedorenko-rois')
    roi_masks_struct = make_fedorenko_roi_masks();
    lr_rois = {roi_masks_struct.left};
else
    lr_rois = cell(1);
end

Nroi = length(lr_rois);
NvoxList = length(voxSizeValues);
Nrepresentation = 2*NvoxList + 1;
Nmodel = 3;
modelLegendName = [string(voxSizeValues), "w2v", string(voxSizeValues)+"+w2v"];
TtestResults_group = zeros(Nrepresentation, Nmodel, 13, Nroi);
TtestResults_allsubj_group = zeros(Nrepresentation, Nmodel, Nroi);
RSAResults_group = zeros(Nrepresentation, Nmodel, 13, Nroi);
for i = 1:length(voxSizeValues)
    load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'RDMAnalysis_results_%s_none_stable%d.mat'), ...
            mode, voxSizeValues(i)));
    TtestResults_group(i, :, :, :) = TtestResults;
    TtestResults_allsubj_group(i, :, :) = TtestResults_allsubj;
    RSAResults_group(i, :, :, :) =RSAResults;
    
    load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'RDMAnalysis_results_%s_word2vec_stable%d.mat'), ...
            mode, voxSizeValues(i)));
    TtestResults_group(NvoxList+1+i, :, :, :) = TtestResults;
    TtestResults_allsubj_group(NvoxList+1+i, :, :) = TtestResults_allsubj;
    RSAResults_group(NvoxList+1+i, :, :, :) =RSAResults;
end
% pure word2vec, voxSizeValues = 0
     load(sprintf(fullfile(BaseDir(), 'group_analysis', ...
            'RDMAnalysis_results_%s_word2vec_stable%d.mat'), ...
            mode, 0));
    TtestResults_group(NvoxList+1, :, :, :) = TtestResults;
    TtestResults_allsubj_group(NvoxList+1, :, :) = TtestResults_allsubj;
    RSAResults_group(NvoxList+1, :, :, :) =RSAResults;
    
%% figure 1: t-test for individual subjects, scatter plot each t value and
% p-value
%{
figure();
hold on;
for mdl = 1:Nrepresentation
        p1 = plot((mdl-0.1) * ones(1, Nsubjects), squeeze(TtestResults_group(mdl, 1, :)), 'o');
        p2 = plot((mdl+0.1) * ones(1, Nsubjects), squeeze(TtestResults_group(mdl, 2, :)), '*');
        p3 = plot((mdl+0.1) * ones(1, Nsubjects), squeeze(TtestResults_group(mdl, 3, :)), 'x');
end
title("T-test on individual subjects' RDM entries");
axis([0, (Nrepresentation+1), 0, 1]);
xticks(1:Nrepresentation);
xticklabels(modelLegendName);
xlabel('number of stable voxels');
ylabel('p-values');
legend([p1, p2, p3], {'taxonomic model', 'thematic model', 'tax + the model'});
saveas(gca, ...
             fullfile(BaseDir(), 'figures', ...
                      'RDMAnalysis_Ttest_eachsubj.png'), 'png');
%}

%% figure 2: t-test for all subjects, plot p-values 
% not that informative right now...
%figure();
%plot(TtestResults_allsubj_group);

%% figure 3: Kendall's correlation
%{
figure();
hold on;
p1 = plot(mean(RSAResults_group(:, 1, :), 3), 'r');
p2 = plot(mean(RSAResults_group(:, 2, :), 3), 'b');
p3 = plot(mean(RSAResults_group(:, 3, :), 3), 'g');

%scatter plot each subject's value
for mdl = 1:Nrepresentation
        plot((mdl-0.15) * ones(1, Nsubjects), squeeze(RSAResults_group(mdl, 1, :)), 'ro');
        plot((mdl) * ones(1, Nsubjects), squeeze(RSAResults_group(mdl, 2, :)), 'bo');
        plot((mdl+0.15) * ones(1, Nsubjects), squeeze(RSAResults_group(mdl, 3, :)), 'go');
end
axis([0, (Nrepresentation+1), -0.05, 0.7]);
xticks(1:Nrepresentation);
xticklabels(modelLegendName);
xlabel('number of stable voxels');
ylabel("Kendall's tau correlation");
legend([p1, p2, p3], {'taxonomic model', 'thematic model',  'tax + the model'});

% do ttest
RSAResults_t = zeros(Nrepresentation, Nmodel);
RSAResults_p = zeros(Nrepresentation, Nmodel);
for i = 1:Nrepresentation
    for mdl = 1:Nmodel
        [h, RSAResults_p(i, mdl), ci, stats] = ttest(RSAResults_group(i, mdl, :), 0, 'Tail', 'right');
        RSAResults_t(i, mdl) = stats.tstat;
    end
end

% draw significant stars
for i = 1:Nrepresentation
    mysigstar(gca, i-0.15, 0.48, RSAResults_p(i, 1));
    mysigstar(gca, i, 0.5, RSAResults_p(i, 2));
    mysigstar(gca, i+0.15, 0.52, RSAResults_p(i, 3));
end
saveas(gca, ...
             fullfile(BaseDir(), 'figures', ...
                      'RDMAnalysis_KendallCorr.png'), 'png');
%}
%% figure 4: bar plot. for each # of stable voxels, plot the correlations groups by ROIs
repre = 3;
for repre = 1:Nrepresentation
    figure('Position', [0, 0, 1000, 700]);
    ax = gca;
    hold on;
    
    % do significance test
    y = squeeze(mean(RSAResults_group(repre, :, :, :), 3))';
    y_std = squeeze(std(RSAResults_group(repre, :, :, :), 0, 3))';
    ts = y ./ (y_std / sqrt(Nsubjects-1));
    ps = tcdf(ts, Nsubjects-1, 'upper');
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
            mysigstar(ax, x(j), ax.YLim(2)-0.01, ps(j, i));
        end
    end
    xticks(1:Nroi);
    roi_names = struct2cell(roi_masks_struct);
    roi_names = squeeze(roi_names(2, :, :));
    xticklabels(roi_names);
    ylabel("Kendall's tau correlation");
    legend(p1, {'taxonomic model', 'thematic model',  'tax + the model'});
    title(sprintf('Correlation for each ROI and each model, at %s voxels\n', modelLegendName{repre}));
    saveas(gca, ...
             sprintf(fullfile(BaseDir(), 'figures', ...
                      'RDMAnalysis_KendallCorr_%s_%s.png'), ...
                      mode, modelLegendName{repre}), 'png');

end
a=1;

