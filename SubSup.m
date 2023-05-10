function [S, H_values, perm_choices, itertime] = SubSup(H, V, perms, submin_algo, maxiter, localmin)
% Compute  min_{S} H(S):= F(S) - G(S) using Submodular-Supermodular procedure of 
% Iyer et al, "Algorithms for Approximate Minimization of the Difference
% Between Submodular Functions, with Applications", 2013.
%
% INPUT
% H: difference of submodular functions H.H(S) = H.F(S)- H.G(S), with
% H.H, H.F and H.G of class SetFct
% perms: choice of heuristic for order of permutation 
%   1- use random order
%   2- order according to decreading gains of G 
%   3- order according to decreading gains of H
%   If more than one perm chosen, try each one at every iter and choose one 
%   that leads to smaller H value
% submin_algo: unconstrained submodular minimization algorithm 
%              returns S, bundle = submin_algo(H, param_H, bundle) = argmin_S H(S)
%              bundle is information needed to warm start algorithm
% maxiter: maximum number of iterations
% localmin: use local minimality as additional stopping criterion

% OUTPUT
% S: argmin H(S)
% Hvalues: H(S) at every iteration

n = length(V);
H_marg_S = zeros(n,1);
G_marg_S = zeros(n,1);
H_values = zeros(maxiter, 1);
itertime = zeros(maxiter, 1);
perm_choices = zeros(maxiter, 1);

S = [];
bundle = cell(length(perms),1);

delete(gcp('nocreate')); % ensure no parallel pool is running
ncpus = str2num(getenv('SLURM_CPUS_ON_NODE'))
if isempty(ncpus) % running on laptop
    parpool(3);
else
    parpool(min(str2num(getenv('SLURM_CPUS_ON_NODE')), 3));
end
% Start the clock.
   tic;
for iter =1:maxiter

    % compute marginals w.r.t to S
    % set current sol of H and G to S, for faster computation of marginals
    if ismember(3, perms)
       [~, H.H] = H.H(S); 
       for i = V
        H_marg_S(i) = add(H.H, S, i) - rmv(H.H, S, i); %H.H(union(S,i), param_H) - H.H(setdiff(S, i), param_H);
       end
    end
    if ismember(2, perms) 
       [~, H.G] = H.G(S);
       for i = V
       G_marg_S(i) = add(H.G, S, i) - rmv(H.G, S, i); %H.G(union(S,i), param_H.G) - H.G(setdiff(S, i), param_H.G);   
       end
    end
    % construct indicator vector of S
    indS = ismember(V, S);
    
    S_next = cell(length(perms),1);
    H_next = zeros(length(perms),1);
    parfor perm = perms
       
        % choose tie breaking order according to greatest gains of G or H w.r.t S
        if perm == 2 
            ties = G_marg_S;
        elseif perm == 3 
            ties = H_marg_S;
        else
            ties = randperm(n);
        end
        subgrad_G = greedy_algo_submodular(indS, H.G, ties);
        
        % choose permutation order according to greatest gains of G or H w.r.t S
%         if perm == 2
%             [~,order] = sort(G_marg_S,'descend');
%         elseif perm == 3 
%             [~,order] = sort(H_marg_S,'descend');
%         else
%             order = V';
%         end
%         
%         % get indices of elements in S in order
%         [~, indS_inorder] = ismember(S, order);
%         % rearrange indices to keep same ordering in order
%         indS_inorder = sort(indS_inorder, 'ascend');
%         % rearrange elements in order to have elements in S first
%         order = [order(indS_inorder); order(setdiff(V, indS_inorder))];
%         %fprintf("Does order contain S after reordering? %d \n", numel(setdiff(order(1:length(S)), S))==0) 
% 
%         % define modular lower bound on G
%         subgrad_G = zeros(n,1);
%         Gval_prev = 0;
%         for i = V
%             Gval = G(order(1:i), param_G);
%             subgrad_G(order(i)) = Gval - Gval_prev; 
%             Gval_prev = Gval;
%         end

        % solve unconstrained sub min problem: min F(A) - subgrad_G(A)
        % G_min_next = @(S,e,AS_pinv) deal(sum(subgrad_G([S;e])), []);
        % S = minimize_WDRsub_projected_subgradient_descent_closure(F,G_min_next,param_H,PGM_param.maxiter,PGM_param.gap_tol, ErrFcns, compute_alphaF,compute_betaG);
        % [H_approx, param_H_approx] = diff_sub_fct(H.F, param_H.F, @(S, param) sum(subgrad_G(S)), struct);
        H_approx = SetFctLinComb({H.F,ModFct(subgrad_G)},[1, -1]);
        if iter == 1
            [S_next{perm}, bundle{perm}] = submin_algo(H_approx);
        else % warm start
            [S_next{perm}, bundle{perm}] = submin_algo(H_approx, bundle{perm});
        end
        H_next(perm) = H.H(S_next{perm});
    end
    
    [H_next, perm_choices(iter)] = min(H_next);
    H_values(iter) = H_next;
    S_next = S_next{perm_choices(iter)};
    
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
end
H_values = H_values(1:iter);
itertime = itertime(1:iter);
perm_choices = perm_choices(1:iter);
end