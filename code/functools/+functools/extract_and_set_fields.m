function [receiver] = extract_and_set_fields(sender, reciever, fields)
% Pull out the named fields from the "sender", and set in "receiver".
% (think of it as a cross-join for struct arrays).
%
% >> x = struct('a', {1 2 3}, 'b', {4 5 6});
% >> y = struct('c', {'one', 'two'});
% >> extract_and_set_fields(x, y, {'a', 'b'})
% 
% ans = 
% 
% 1x6 struct array with fields:
% 
%     c
%     a
%     b


if ~isstruct(sender), error('Sender must be struct array'); end
if ~isstruct(reciever), error('Receiver must be struct array'); end

if ~iscellstr(fields)
    if ischar(fields)
        fields = {fields};
    else
        error('Fields must be cellstr');
    end
end

import functools.*;

receiver = map(@(send) ...
               functools.map(make_extractor_setter(send, fields), ...
                             reciever), ...
               sender);
receiver = cell2mat([receiver{:}]);

end

function fcn = make_extractor_setter(send, fields)
funs = functools.map(@(field) @(rec) setfield(rec, field, send.(field)), fields);
fcn = functools.compose(funs{:});
end
