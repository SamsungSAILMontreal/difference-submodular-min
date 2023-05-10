function [counts] = data_count(data,n_rows)
%Count the observed frequency of permutations in the data set (data).
%For example, if "[0, 0, 1]" is the value contained in some row of data,
%then we want to count how many times this permutation is observed in data.
    [data_uniq,ia,ic] = unique(data, 'rows', 'stable');
    counts = accumarray(ic, 1);

end

