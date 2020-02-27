function RDMDifferenceAnalysis(sessions, mode, wordModel, NstableVox)
%1. gather beta matrices and RDMs (loadBetaAndRDM function)
% 2. create RDM based on fMRI subtraction of word pairs
% 3. plot them and see what happens

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

% 2. create the category difference RDM
% for each theme, take the pattern difference between tool and person words
% do the same for other two

Nsubj = length(RDM);
betaDiff = cell(size(beta));
categoryDiffRDM = cell(size(beta));
sumCategoryDiffRDM = zeros(45);
% implicitly assumes stimEntries were ordered according to category and
% theme
categoryComb = (nchoosek(1:3, 2) - 1) * 15;

for subj = 1:Nsubj
    ctr = 1;
    for cate_idx = 1:size(categoryComb, 1)
        for theme = 1:15
            betaDiff{subj}(ctr, :) = beta{subj}(theme + categoryComb(cate_idx, 1), :) - ...
                                                 beta{subj}(theme + categoryComb(cate_idx, 2), :);
            ctr = ctr + 1;
        end
    end
    % calculate RDM
    categoryDiffRDM{subj} = corr(betaDiff{subj}');
    sumCategoryDiffRDM = sumCategoryDiffRDM + categoryDiffRDM{subj};
end

% form entries within each pattern difference pair
triu_idx = find(triu(ones(15), 1) == 1);
withinBetweenGroup = cell(Nsubj, 2);
for subj = 1:Nsubj
    for gp1 = 1:3
        for gp2 = (gp1):3
            entries = categoryDiffRDM{subj}( ...
                    (gp1-1)*15 + (1:15), (gp2-1)*15 + (1:15));
            if gp1 == gp2 % same group
                withinBetweenGroup{subj, 1} = [ ...
                    withinBetweenGroup{subj, 1}; entries(triu_idx)];
            else
                withinBetweenGroup{subj, 2} = [ ...
                    withinBetweenGroup{subj, 2}; entries(triu_idx)];
            end
        end
    end
end
 
% do t-test on each subject
p_values = zeros(Nsubj, 1);
for subj = 1:Nsubj
    [h, p_values(subj)] = ttest(withinBetweenGroup{subj, 1}, ...
                               withinBetweenGroup{subj, 2}, 'tail', 'Right');
end

% average for each subject, and do group-level test
avgGroup = zeros(Nsubj, 2);
for subj = 1:Nsubj
    avgGroup(subj, 1) = mean(withinBetweenGroup{subj, 1});
    avgGroup(subj, 2) = mean(withinBetweenGroup{subj, 2});
end
[h, p_value_group] = ttest(avgGroup(:, 1) - avgGroup(:, 2), 0, 'tail', 'Right');
a=1;
