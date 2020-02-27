function [c] = cellify(x)
% turn anything into a cell array, with the power of map.

c = functools.map(@(xx) xx, x);