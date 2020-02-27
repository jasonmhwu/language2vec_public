function makeConditionsFiles(subjectNum, prepost, output_dir_template, matfile_template)
%%% makeConditionsFiles(subjectNum, prePost)
%%%
%%% Creates condition timing files for SPM analysis, all trials in a single
%%% regressor so we can compare trials to baseline 

fprintf('Generating conditions files for subject %d, %s-scan\n', ...
    subjectNum, prepost);

% get runs and write some info to the console
runs = parseOneSession(subjectNum, prepost);
nruns = length(runs);

fprintf('Found data for %d runs:\n', nruns);
for i = 1:nruns
    runInfo = cellfun(@(name) getColumn(getRows(runs{i}, 1), name), ...
        {'task', 'Time', 'Date'});
    runInfo{4} = length(getColumn(runs{i}, 'Date'));
    fprintf('  Run %2d: %15s at %s, %s (%d trials)\n', i, runInfo{:});
end

% get trials per run and total number for subject
trialsPerRun = cellfun(@(run) length(run.data{1}), runs);
totalTrials = sum(trialsPerRun);

% initialize conditions struct for each run separately (to avoide empty
% conditions):  
%   onsets is a runs x conditions cell array
%   durations is a runs x conditions cell array (of zeros)
%   names is a runs x conditions cell array of condition names
trialConditions = struct();
for runNum = 1:nruns
    run = runs{runNum};
    trials = trialsPerRun(runNum);
    trialConditions(runNum).onsets = num2cell(getColumn(run, 'trialSecs'));
    trialConditions(runNum).durations = num2cell(zeros(trials, 1));
    
    nameData = [getColumn(run, 'task'), ...
                num2cell(getColumn(run, 'block')), ...
                repmat({'trial'}, [trials, 1]), ...
                num2cell(getColumn(run, 'trial'))];
            
    trialConditions(runNum).names = cell(trials, 1);
    for trial = 1:trials
        trialConditions(runNum).names{trial} = ...
            sprintf('%s-%d-%s-%d', nameData{trial, :});
    end
end



%% create mat files
% default output dir is: 
% ...imaging-analysis/models/<model dir>/subjects/Lin00n/<prepost>/conditions
if nargin < 3
    output_dir_template = fullfile(animalImagingBaseDir(), ...
                                   modelSubdir(), ...
                                   'subjects', ...
                                   'Lin%03d', ...
                                   '%s', ...
                                   'conditions');
end
output_dir = sprintf(output_dir_template, subjectNum, prepost);

% default file name template is 
if nargin < 4
    matfile_template = 'conditions_run%03d.mat';
end

fprintf('Writing condition .mat files in %s\n', output_dir);
                  
if ~ exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:nruns
    
    onsets = trialConditions(i).onsets;
    durations = trialConditions(i).durations;
    names = trialConditions(i).names;

    save(fullfile(output_dir, sprintf(matfile_template, i)), ...
        'onsets', 'durations', 'names');
end