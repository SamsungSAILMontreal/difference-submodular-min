function [local_minimality, U] = check_localmin(H, S, V)
% Check if set S is a local minimum of H
% return U=S if S is a local min, otherwise return U with smallest H(U) 
% s.t U differs from S by just one element
% local_minimality = H(S) - min{H(U) : U differs from S by one element}

% set current sol of H to solution, for faster computation of marginals
[H_S, H.H] = H.H(S);
min_H_local = inf;
for i= setdiff(V, S)
    H_Si = add(H.H, S , i);
    if H_Si < min_H_local
        min_H_local = H_Si;
        U = union(S,i);
    end
end
for i= S(:)' 
    H_Smi = rmv(H.H, S , i);
    if H_Smi < min_H_local
        min_H_local = H_Smi;
        U = setdiff(S,i);
    end
end
local_minimality = H_S - min_H_local;
if local_minimality <= 0
    U = S;
end
end