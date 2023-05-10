function [S, H_values, itertime] = ModMod(H, V, modFcts, perms, maxiter, localmin)
% Compute  min_{S} H(S):= F(S) - G(S) using Modular-Modular procedure of 
% Iyer et al, "Algorithms for Approximate Minimization of the Difference
% Between Submodular Functions, with Applications", 2013.
%
% INPUT
% H: difference of submodular functions H.H(S) = H.F(S)- H.G(S), with
% H.H, H.F and H.G of class SetFct
% V: ground set
% modFct: choice of modular upper bound on F 
%   1- use first modular upper bound in Iyer el al, 2013
%   2- use second modular upper bound in Iyer el al, 2013
%   If more than one modFct chosen, try each one at every iter and choose 
%   one that leads to smaller H value
% perm: choice of heuristic for order of permutation 
%   1- use random order
%   2- order according to decreading gains of G 
%   3- order according to decreading gains of H
%   If more than one perm chosen, try each one at every iter and choose one 
%   that leads to smaller H value
% maxiter: maximum number of iterations
% localmin: use local minimality as additional stopping criterion
% OUTPUT
% S: argmin H(S)
% Hvalues: H(S) at every iteration
% itertime: elapsed time after each iteration
n = length(V);
Fsingletons = zeros(n,1);
FVsingletons = zeros(n,1);
F_marg_S = zeros(n,1);
G_marg_S = zeros(n,1);
H_values = zeros(maxiter, 1);
itertime = zeros(maxiter, 1);
[FV_val, FV_obj] = H.F(V);
[F0_val, F0_obj] = H.F([]); % should be zero

for i = V
    if ismember(1, modFcts)
       Fsingletons(i) = add(F0_obj,[],i) - F0_val;
    end
    if ismember(2, modFcts)
       FVsingletons(i) = FV_val - rmv(FV_obj, V, i);
    end
end

S = [];
Sc = V;

delete(gcp('nocreate')); % ensure no parallel pool is running
ncpus = str2num(getenv('SLURM_CPUS_ON_NODE'))
if isempty(ncpus) % running on laptop
    parpool(4);
else
    parpool(min(str2num(getenv('SLURM_CPUS_ON_NODE')), 6));
end
% Start the clock.
   tic;
for iter =1:maxiter
       
    % compute marginals w.r.t to S
    % set current sol of F and G to S, for faster computation of marginals
    [~, H.F] = H.F(S); 
    [~, H.G] = H.G(S);
    for i = V
       F_marg_S(i) = add(H.F, S, i) - rmv(H.F, S, i); %F(union(S,i)) - F(setdiff(S, i));
       G_marg_S(i) = add(H.G, S, i) - rmv(H.G, S, i); %G(union(S,i)) - G(setdiff(S, i));       
    end
    % construct indicator vector of S
    indS = ismember(V, S);
    
    config_size = [length(perms),length(modFcts)];
    nconfigs = prod(config_size);
    S_next = cell(nconfigs,1);
    H_next = zeros(nconfigs,1);
    parfor config = 1:nconfigs
        [perm, modFct] = ind2sub(config_size, config);
        
        % define modular upper bound on F
        M_F = zeros(n,1);
        if modFct == 1
           M_F(S) = - F_marg_S(S);
           M_F(Sc) = Fsingletons(Sc);
        else
           M_F(S) = - FVsingletons(S);
           M_F(Sc) = F_marg_S(Sc);
        end

         % choose tie breaking order according to greatest gains of G or H w.r.t S
        if perm == 2 
            ties = G_marg_S;
        elseif perm == 3 
            ties = F_marg_S - G_marg_S;
        else
            ties = randperm(n);
        end
        subgrad_G = greedy_algo_submodular(indS, H.G, ties);

        S_next{config} = find(M_F - subgrad_G < 0);
        H_next(config) = H.H(S_next{config});
    end
    
    [H_next, ind] = min(H_next);
    H_values(iter) = H_next;
    S_next = S_next{ind};
    
    itertime(iter) = toc;
    if isempty(setxor(S, S_next))
        if localmin
            [local_minimality, U] = check_localmin(H, S, V);
            if local_minimality>0
                fprintf("Algorithm converged after %d iterations but not to a local min, local minimality = %e \n", iter, local_minimality)
                S_next = U;
                H_values(iter) = H_values(iter) - local_minimality;
            else
                fprintf("Algorithm converged after %d iterations to a local min, local minimality = %e \n", iter, local_minimality)
                break
            end
        else
            fprintf("Algorithm converged after %d iterations \n", iter)
            break
        end
    end
    S = S_next;
    Sc = setdiff(V, S);
end
H_values = H_values(1:iter);
itertime = itertime(1:iter);
end