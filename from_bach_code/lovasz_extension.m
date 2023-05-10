function f = lovasz_extension(w,F)
% compute lovasz extension
f = w'* greedy_algo_submodular(w,F);