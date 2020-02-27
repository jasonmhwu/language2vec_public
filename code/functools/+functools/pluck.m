function [out] = pluck(in, varargin)
% pull out nested struct fields from struct arrays.  god help you if the
% nested fields can't be concatenated.
% 
% Returns a flattened array:
%
% >> x = struct('a', ..
% >>            { struct('b', {1, 2, 3} ) , ...
% >>              struct('b', {4, 5} ) } );
% >> pluck(x, 'a', 'b')
% 
% ans = 
%      1    2    3    4    5
% 
% Why does this exist at all?  Because you might want to do something like 
% 
% >> x.a.b
% 
% or
% 
% >> [x.a].b
% 
% but you're not allowed to. :(

if nargin < 1
    error('Gotta have something to pluck');
end

if ~ isstruct(in)
    error('Can only pluck from a struct array');
end

top_field = varargin{1};
rest_fields = varargin(2:end);

if isempty(rest_fields)
    % no other levels remaining? return cell array of field entries
    out = { in.(top_field) }';
else
    % recurse
    out = functools.pluck([in.(top_field)], rest_fields{:});
end
