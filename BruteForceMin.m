function [S_min, H_min] = BruteForceMin(H, V)
 
n = length(V);
S_min = [];
H_min = H(S_min);

% Go over all possible subsets
for i=1:2^n-1
    S = dec2bin(i,n) - '0';
    S = V(S~=0);
    H_S = H(S);
    if H_S < H_min
        H_min = H_S;
        S_min = S;
    end
end
end