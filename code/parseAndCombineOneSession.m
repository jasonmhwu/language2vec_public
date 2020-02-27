function all_runs = parseAndCombineOneSession(subjectNumber, varargin)
% combine multiple runs of behavioral data from parseOneSession

% read in individual runs.  each run is read in as a "data frame" type
% structure: each row is an observation, each column is a variable (date,
% time, animal_x, etc.) and each row is an observation (trial).
runs = parseOneSession(subjectNumber);

% runs is a cell array of each individual runs, and they need to be assembled
% into one big data frame using `addRow`.
all_runs = addColumn(runs{1}, runs{1}.runs, 'runNum');
for run_num = 2:length(runs)
    all_runs = addRow(all_runs, addColumn(runs{run_num}, runs{run_num}.runs, 'runNum'));
end

all_runs = rmfield(all_runs, 'runs');