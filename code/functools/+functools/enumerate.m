% like map, but first argument to fn is index, and second is value.

function retval = enumerate(fn, args)

args_cell = functools.cellify(args);

idxs_in_shape = zeros(size(args));
idxs_in_shape(:) = 1:numel(args);

retval = functools.map(@(idx) fn(idx, args_cell{idx}), idxs_in_shape);
