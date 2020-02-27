function visualizeMDS( stimMatrix, colName, figname)
%calculate RDM for the specified column, and color-code them based on
%type (colr) and semantic category (textLabel)
if nargin < 3
    figname = '';
end
    

distance = 1 - corr(stimMatrix.(colName)');
mds_vector = cmdscale(distance, 2);

figure('Position', [10, 10, 1000, 1000]);
hold on;

uniqColor = unique(stimMatrix.wordType);
colorMap = {'r', 'g', 'b'};
for col_idx = 1:length(uniqColor)
    idx = strcmp(stimMatrix.wordType, uniqColor{col_idx});
    text(mds_vector(idx, 1), mds_vector(idx, 2), stimMatrix.Stimulus(idx), ...
        'Color', colorMap{col_idx}, 'FontSize', 12);
end
axis([min(mds_vector(:, 1))-0.15, max(mds_vector(:, 1))+0.15, min(mds_vector(:, 2))-0.15, max(mds_vector(:, 2))+0.15]);
hold off;

% save figure
if ~isempty(figname)
    saveas(gcf, figname);
end
close();
end
