function [S, H_values, itertime] = SupSub(H, V, modFcts, submax_algo, maxiter, localmin)
% Compute  min_{S} H(S):= F(S) - G(S) using Supermodular-Submodular procedure of 
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
% submax_algo: unconstrained submodular minimization algorithm 
% returns S = submax_algo(H, H_marg) = approximation of argmax_S H(S)
% maxiter: maximum number of iterations
% localmin: use local minimality as additional stopping criterion

% OUTPUT
% S: argmin H(S)
% Hvalues: H(S) at every iteration

n = length(V);
Fsingletons = zeros(n,1);
FVsingletons = zeros(n,1);
F_marg_S = zeros(n,1);
H_values = zeros(maxiter, 1);
itertime = zeros(maxiter, 1);

[FV_val, FV_obj] = H.F(V);
[F0_val, F0_obj] = H.F([]);

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
    parpool(2);
else
    parpool(min(str2num(getenv('SLURM_CPUS_ON_NODE')), 2));
end
% Start the clock.
   tic;
for iter =1:maxiter

    % compute marginals w.r.t to S
    % set current sol of F to S, for faster computation of marginals
    [~, H.F] = H.F(S); 
    for i = V
       F_marg_S(i) = add(H.F, S, i) - rmv(H.F, S, i); %F(union(S,i)) - F(setdiff(S, i));
    end
    
    S_next = cell(length(modFcts),1);
    H_next = zeros(length(modFcts),1);
    
    parfor modFct = modFcts
        % define modular upper bound on F
        M_F = zeros(n,1);
        if modFct == 1
           M_F(S) = - F_marg_S(S);
           M_F(Sc) = Fsingletons(Sc);
        else
           M_F(S) = - FVsingletons(S);
           M_F(Sc) = F_marg_S(Sc);
        end
        
        % solve unconstrained sub max problem: max G(A) - M_F(A)
        % Hmax = @(A) H.G(A, param_H.G) - sum(M_F(A));
        % Hmax_marg = @(A,i) H.G([A,i], param_H.G)- H.G(A, param_H.G) - M_F(i);
        Hmax = SetFctLinComb({H.G,ModFct(M_F)},[1, -1]);
        S_next{modFct} = submax_algo(Hmax);
        H_next(modFct) = H.H(S_next{modFct});
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