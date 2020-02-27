function [classes] = make_classes_wordtype(session)
% generate classes for decoding: person vs. build vs. tool
% 
% input:
%   session: info struct for one session,
%   classes: row: label, column: identity for cross validation selection
    
    all_runs = parseAndCombineOneSession(session.SubjectNum);
    actualTrials = getColumn(all_runs, 'actualTrials');
    for i = 1:length(all_runs.varnames)
        all_runs.data{i} = all_runs.data{i}(actualTrials);
    end
    wordType = getColumn(all_runs, 'wordType');
    stimCategory = getColumn(all_runs, 'stimCategory');
    uniqueVal = {'person', 'building', 'tool'};
    classes = cell(3, 15);
    for type = 1:length(uniqueVal)
        for stimCat = 1:15
        classes{type, stimCat} = find(...
            strcmp( wordType, uniqueVal(type)) & stimCategory == stimCat);
        end
    end
end