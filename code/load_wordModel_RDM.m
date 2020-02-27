function [RDM, stimMatrix] = load_wordModel_RDM(stimMatrix, model)
% load word2vec or glove RDM representation from mat file
if nargin < 2
    model = 'word2vec';
end

if strcmp(model, 'word2vec')
    mat = matfile('/home/URMC-SH/mwu34/Desktop/word2vec_project/word2vec_representation.mat');
    whosmat = whos(mat);
else    
end

for i = 1:size(stimMatrix, 1)
    try
        stimMatrix.wordModel_vec{i} = mat.(stimMatrix.Stimulus{i});
    catch
        if ismember(stimMatrix.Stimulus{i}(1), 'a':'z')
            stimMatrix.wordModel_vec{i} = mat.([upper(stimMatrix.Stimulus{i}(1)), lower(stimMatrix.Stimulus{i}(2:end))]);
        else
            stimMatrix.wordModel_vec{i} = mat.(lower(stimMatrix.Stimulus{i}));
        end
    end
end
stimMatrix.wordModel_vec = double(cell2mat(stimMatrix.wordModel_vec));
RDM = corr(stimMatrix.wordModel_vec');
