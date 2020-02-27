function [grouped, group_levels] = group_by(group_fn, x)

group_var = functools.reduce(@vertcat, functools.map(group_fn, x));
group_levels = unique(group_var);
try 
    grouped = functools.map(@(level) x(strcmp(level, group_var)), group_levels);
catch
    grouped = functools.map(@(level) x(group_levels == level), group_levels);
end
