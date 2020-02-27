function [split_struct, split_labels] = split(input, split_on)
% Split struct array input on (one or more) named fields. If split_on names
% more than one field, struct will be split based on all unique combinations of
% values of those fields
% 
% [split_struct, split_labels] = split(input, split_on)
% 
% Input:
%   input: struct array to split
%   split_on: cell array of field names to split on.
% Output:
%   split_struct: cell array of splits. Each element is subset of input sharing
%      one combination of values of fields in `split_on`.
%   split_labels: struct array given values of split_on for each element of
%     `split_struct`.

import functools.*;

% construct an index variable by pasting fields together
paste_fields = @(fields) functools.reduce(@(x,y) [x '_' y], fields);
stringify_field = @(field) if_(isnumeric(field), ...
                               @() num2str(field), ...
                               @() field);
stringify_fields = @(fields) functools.map(stringify_field, fields);
extract_fields = @(row) functools.map(@(field) row.(field), split_on);

index_ = map(functools.compose(extract_fields, ...
                               stringify_fields, ...
                               paste_fields), ...
             input);

split_struct = map(@(i) input(strcmp(i, index_)), unique(index_));

split_labels = map(@(one) functools.extract_and_set_fields(one(1), ...
                                                           struct(), ...
                                                           split_on), ...
                   split_struct);
split_labels = cell2mat(split_labels);
