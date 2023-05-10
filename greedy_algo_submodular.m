 function [x,Fvalues,I] = greedy_algo_submodular(w,F, ties)
% Given a submodular function F of class SetFct 
% returns x=gradient of the lovasz extension
% x \in argmax w^Ts , s \in B(F)
% Fvalues=values of F(j1 .. jk), where w_j1>= ...>= w_jk
% I= j1 ... jn (sorted coefficients of w)
% ties: values to use to break ties when sorting w

n = length(w);
x = zeros(n,1); 
if nargin<3
    [~,I] = sort(w,'descend');
else
    [~,I] = sortrows([w(:), ties(:)], 'descend');
end
    
Fvalues = zeros(n,1);

[Fvalues(1), F] = F(I(1));
x(I(1)) = Fvalues(1);
for i=2:n
    [Fvalues(i), F] = add(F, I(1:i-1), I(i)); 
    x(I(i)) = Fvalues(i) - Fvalues(i-1);
end

