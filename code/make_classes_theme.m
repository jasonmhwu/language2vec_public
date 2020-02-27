function [classes] = make_classes_theme(session)
% generate classes for decoding: person vs. build vs. tool
% 
% input:
%   session: info struct for one session,
%   classes: { [person trials], [building trials], [tool trials] }
    
    all_runs = parseAndCombineOneSession(session.SubjectNum);
    actualTrials = getColumn(all_runs, 'actualTrials');
    for i = 1:length(all_runs.varnames)
        all_runs.data{i} = all_runs.data{i}(actualTrials);
    end
    wordType = getColumn(all_runs, 'wordType');
    stimCategory = getColumn(all_runs, 'stimCategory');
    uniqueVal = {'person', 'building', 'tool'};
    classes = cell(15, 3);
    for type = 1:length(uniqueVal)
        for stimCat = 1:15
        classes{stimCat, type} = find(...
            strcmp( wordType, uniqueVal(type)) & stimCategory == stimCat);
        end
    end
end