% check for head motion artifacts
plotTimeSeries = 1;

subjects = subj_info();
for subj_idx = 1:length(subjects)
    subj = subjects(subj_idx);
    folderName = fullfile(BaseDir(), 'data-imaging/', ...
        subj.ScanningDir, subj.ScanningDate, ...
        'Mcverter_Dicom_conversion');
    cd(folderName);
    textName = dir('*/rp*.txt');
    if plotTimeSeries
        figure('Position', [0, 0, 700, 700]);
        title(subj.ScanningDir);
    end
    for i = 1:length(textName)
        fid = fopen(fullfile( ...
            textName(i).folder, textName(i).name), 'r');
        sizeArray = [6 Inf];
        timeSeries = fscanf(fid, '%f', sizeArray);
        if plotTimeSeries
            subplot(length(textName), 1, i);
            plot(timeSeries(1:3, :)');
        end
    end
end