function col = getColumn(data, colname)

% get indices of matching column names
whichCol = strcmp(colname, data.varnames);

% throw error if not found
if ~any(whichCol) 
    error('Tried to access non-existent column %s', colname);
end

% return requested column(s)
if sum(whichCol)==1
    col = data.data{whichCol};
else
    col = data.data(whichCol);
end

end