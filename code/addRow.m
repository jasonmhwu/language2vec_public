function data = addRow(data, rows)
% function data = addRow(data, rows)
%
% Add a row (or rows) to an existing data frame.  New rows need to be
% provided as a data frame struct, with fields data and varnames.

if isfield(rows, 'varnames') && isfield(rows, 'data')
    % input is a data frame itself; iterate over varnames, looking for
    % matching one in existing data and adding new rows
    for varname = rows.varnames
        var_index_old = find(strcmp(varname, data.varnames));
        var_index_new = find(strcmp(varname, rows.varnames));
        if var_index_old
            data.data{var_index_old} = [data.data{var_index_old}; rows.data{var_index_new}];
        else
            error('Variable not found: %s', varname);
        end
    end
else
    error('Need to provide rows in data frame format (struct with data and varnames fields)');
end