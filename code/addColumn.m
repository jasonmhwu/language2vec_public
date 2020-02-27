function data = addColumn(data, column, colname)
% function data = addColumn(data, column, colname)
%
% adds a column to an existing data structure, `data`.  `column` can either
% be a vector (or cell vector) or a scalar (which will be broadcasted into
% a vector of the same size as the first column already in `data`).
% returns the updated data structure.

% broadcast singletons
if length(column)==1
    column = repmat(column, size(data.data{1}));
end

data.data{1, end+1} = column(:);
data.varnames{1, end+1} = colname;

end