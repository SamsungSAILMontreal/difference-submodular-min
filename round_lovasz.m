function [S, F_S] = round_lovasz(x, F, Fvalues, order)
% round fractional solution x to set S s.t. F(S)<= f_L(x)
% optional args:
% Fvalues: values of F(j1 .. jk), where x_j1>= ...>= x_jk
% order: j1, ..., jn
if nargin < 4
    [~,Fvalues,order] = greedy_algo_submodular(x,F);
end
% empty set is also a sublevel
[F_S,ind_min] = min([Fvalues; 0]);
if ind_min > length(Fvalues)
    S = [];
else
    S = order(1:ind_min);
end
end