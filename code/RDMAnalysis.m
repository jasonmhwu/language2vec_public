function RDMAnalysis(sessions, mode, wordModel, NstableVox)
%1. gather beta matrices and RDMs (loadBetaAndRDM function)
% 2ab. t-test: group entries according to taxonomic or thematic categories
% 2c. RSA: create model RDM for taxonomic or thematic categories

if nargin < 1
    sessions = subj_info();
end
if nargin < 2
    mode = 'whole-brain';
end
if nargin < 3
    wordModel = 'word2vec'; % 'none' means not appending word2vec
end
if nargin < 4
    NstableVox = 100;
end

output_dir = fullfile(BaseDir(), 'group_analysis');

% 1. gather beta matrices and RDMs
[beta, RDM, stimEntries] = loadBetaAndRDM(sessions, mode, wordModel, NstableVox);

% 2a. sort entries according to taxonomic or thematic categories
thematicCorr = cell(size(RDM));
taxonomicCorr = cell(size(RDM));
thematicNonCorr = cell(size(RDM));
taxonomicNonCorr = cell(size(RDM));
sumTaxTheCorr = cell(size(RDM));
sumTaxTheNonCorr = cell(size(RDM));

taxonomicRDM = zeros(size(RDM{1}));
thematicRDM = zeros(size(RDM{1}));
sumTaxTheRDM = zeros(size(RDM{1}));
for i = 1:size(stimEntries, 1)
    for j = (i+1) : size(stimEntries, 1)
        if strcmp(stimEntries.wordType{i}, stimEntries.wordType{j})
            taxonomicRDM(i, j) =1;
            sumTaxTheRDM(i, j) = 1;
        else
            taxonomicRDM(i, j) = -1;
            sumTaxTheRDM(i, j) = -1;
        end
        if stimEntries.stimCategory(i) == stimEntries.stimCategory(j)
            thematicRDM(i, j) = 1;
            sumTaxTheRDM(i, j) = 1;
        else
            thematicRDM(i, j) = -1;
        end
    end
end
for sess_idx = 1:size(RDM, 1)
    for roi_idx = 1:size(RDM, 2)
        taxonomicCorr{sess_idx, roi_idx} = ...
            RDM{sess_idx, roi_idx}(taxonomicRDM == 1);
        taxonomicNonCorr{sess_idx, roi_idx} = ...
            RDM{sess_idx, roi_idx}(taxonomicRDM == -1);
        thematicCorr{sess_idx, roi_idx} = ...
            RDM{sess_idx, roi_idx}(thematicRDM == 1);
        thematicNonCorr{sess_idx, roi_idx} = ...
            RDM{sess_idx, roi_idx}(thematicRDM == -1);
        sumTaxTheCorr{sess_idx, roi_idx} = ...
            RDM{sess_idx, roi_idx}(sumTaxTheRDM  == 1);
        sumTaxTheNonCorr{sess_idx, roi_idx} = ...
            RDM{sess_idx, roi_idx}(sumTaxTheRDM == -1);
    end
end

% 2b1. for each subject, do a t-test on both taxonomic and thematic groups
TtestResults = zeros(3, size(RDM, 1), size(RDM, 2));
for subj = 1:length(sessions)
    for roi = 1:size(RDM, 2)
        [h, TtestResults(1, subj, roi)] = ttest2(taxonomicCorr{subj, roi}, taxonomicNonCorr{subj, roi}, 'Tail', 'right');
        [h, TtestResults(2, subj, roi)] = ttest2(thematicCorr{subj, roi}, thematicNonCorr{subj, roi}, 'Tail', 'right');
        [h, TtestResults(3, subj, roi)] = ttest2(sumTaxTheCorr{subj, roi}, sumTaxTheNonCorr{subj, roi}, 'Tail', 'right');
    end
end

%2b2. group coefficients from all subjects together
TtestResults_allsubj = zeros(3, size(RDM, 2));
for roi = 1:size(RDM, 2)
        [h, TtestResults_allsubj(1, roi)] = ttest2(reshape(cell2mat(taxonomicCorr(:, roi)), [], 1), ...
            reshape(cell2mat(taxonomicNonCorr(:, roi)), [], 1), 'Tail', 'right');
        [h, TtestResults_allsubj(2, roi)] = ttest2(reshape(cell2mat(thematicCorr(:, roi)), [], 1), ...
            reshape(cell2mat(thematicNonCorr(:, roi)), [], 1), 'Tail', 'right');
        [h, TtestResults_allsubj(3, roi)] = ttest2(reshape(cell2mat(sumTaxTheCorr(:, roi)), [], 1), ...
            reshape(cell2mat(sumTaxTheNonCorr(:, roi)), [], 1), 'Tail', 'right');
end

% 2c. for each subject, make a binary vector that has 0 with noncorr group
% and 1 with corr group, do Kendall's tau correlation, and do final t-test
RSAResults = zeros(3, size(RDM, 1), size(RDM, 2));
for subj = 1:length(sessions)
    for roi = 1:size(RDM, 2)
        tax_roiRDMVec = [taxonomicCorr{subj, roi}', taxonomicNonCorr{subj, roi}'];
        tax_modelRDMVec = [ones(1, length(taxonomicCorr{subj, roi})), zeros(1,  length(taxonomicNonCorr{subj, roi}))];
        RSAResults(1, subj, roi) = corr(tax_roiRDMVec', tax_modelRDMVec', 'Type', 'Kendall');
        
        the_roiRDMVec = [thematicCorr{subj, roi}', thematicNonCorr{subj, roi}'];
        the_modelRDMVec = [ones(1, length(thematicCorr{subj, roi})), zeros(1,  length(thematicNonCorr{subj, roi}))];
        RSAResults(2, subj, roi) = corr(the_roiRDMVec', the_modelRDMVec', 'Type', 'Kendall');
        
        sum_roiRDMVec = [sumTaxTheCorr{subj, roi}', sumTaxTheNonCorr{subj, roi}'];
        sum_modelRDMVec = [ones(1, length(sumTaxTheCorr{subj, roi})), zeros(1,  length(sumTaxTheNonCorr{subj, roi}))];
        RSAResults(3, subj, roi) = corr(sum_roiRDMVec', sum_modelRDMVec', 'Type', 'Kendall');
    end
end
if size(RDM, 2) == 1 % whole-brain, report results here
    [h, p, ci, stats] = ttest(RSAResults(1, :), 0, 'Tail', 'right');
    fprintf('taxonomic model RSA group t-test: mean Kendall correlation = %f, t=%f, p=%f\n', ...
        mean(RSAResults(1, :)), stats.tstat, p);
     [h, p, ci, stats] = ttest(RSAResults(2, :), 0, 'Tail', 'right');
    fprintf('thematic model RSA group t-test: mean Kendall correlation = %f, t=%f, p=%f\n', ...
        mean(RSAResults(2, :)), stats.tstat, p);
    [h, p, ci, stats] = ttest(RSAResults(3, :), 0, 'Tail', 'right');
    fprintf('sum taxonomic and thematic model RSA group t-test: mean Kendall correlation = %f, t=%f, p=%f\n', ...
        mean(RSAResults(3, :)), stats.tstat, p);
else % harvard-oxford, do t-test on each ROI, and output uncorrected significant resutls
    RDM_t = zeros(3, size(RDM, 2));
    RDM_p = zeros(3, size(RDM, 2));
    for roi = 1:size(RDM, 2)
        for mdl = 1:3
             [h, RDM_p(mdl, roi), ci, stats] = ttest(RSAResults(mdl, :, roi), 0, 'Tail', 'right');
             RDM_t(mdl, roi) = stats.tstat;
        end
    end
    % test significant ROIs here
    % example: find(RDM_p(1, :) < 0.05)
end

% finally, store results in cache
save(sprintf(fullfile(output_dir, 'RDMAnalysis_results_%s_%s_stable%d.mat'), mode, wordModel, NstableVox), ...
    'taxonomicCorr', 'taxonomicNonCorr', 'thematicCorr', 'thematicNonCorr', ...
    'TtestResults', 'TtestResults_allsubj', 'RSAResults');

    
a=1;
