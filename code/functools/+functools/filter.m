function out = filter(predicate, in)
% Applies a function to an input array, and returns the sub-array where
% predicate evalutes to truthy.
% 
% Example:
% 
% >> is_odd = @(x) mod(x,2) != 0;
% >> odd_numbers = filter(is_odd, 1:10);

matches_predicate = logical(cell2mat(functools.map(predicate, in)));
out = in(matches_predicate);
