function [S, x, H_values, hL_values, gaps, perm_choices, itertime] = regDC(H, V, rho, perms, prox, maxiter, tol, warm_start, round, acc, q, localmin)
% Compute  min_{S} H(S):= F(S) - G(S)
% by applying DC algo to the continuous problem
% min_x f_L(x) + delta_[0,1]^n(x) + rho ||x||^2 - g_L(x) - rho ||x||^2
% where f_L and g_L are the Lovasz extensions of F and G
% INPUT
% H: difference of submodular functions H.H(S) = H.F(S)- H.G(S), with
% H.H, H.F and H.G of class SetFct
% V: ground set
% rho: regularization parameter
% perm: choice of heuristic for order to use to break ties when computing
% subgradient of g_L
%   1- use random order
%   2- order according to decreading gains of G 
%   3- order according to decreading gains of H
%   If more than one perm chosen, try each one at every iter and choose one 
%   that leads to smaller H value
% prox: proximal operator of f_L restricted to [0,1]^n returns
% w = prox(F,param_F,y,rho, w_init) =  argmin_{w in [0, 1]^n} f_L(w) + rho/2 ||w||^2 - w^Ty
% submin_algo: submodular minimization algorithm returns S = submin_algo(H, param_H)
% maxiter: maximum number of iterations
% localmin: use local minimality as additional stopping criterion

% OUTPUT
% S: argmin H(S)
% Hvalues: H(S) at every iteration

n = length(V);
x = zeros(n,1);
S = [];
z = x;
x_prev = x;
t = 1;

H_values = zeros(maxiter, 1);
hL_values = zeros(maxiter, 1);
gaps = zeros(maxiter, 1);
perm_choices = zeros(maxiter, 1);
itertime = zeros(maxiter, 1);
H_marg_S = zeros(n,1);
G_marg_S = zeros(n,1);

delete(gcp('nocreate')); % ensure no parallel pool is running
ncpus = str2num(getenv('SLURM_CPUS_ON_NODE'))
if isempty(ncpus) % running on laptop
    parpool(3);
else
    parpool(min(ncpus, 3));
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
    
    if acc
        t_next = (1 + sqrt(1 + 4*t^2))/2;
        z = x + (t-1)*(x - x_prev)/t_next;
        t = t_next;
    end

    if acc && lovasz_extension(z, H.H) <= max(hL_values(max(iter-q,1):iter))
        v = z;
    else
        v = x;
    end
    
    x_next = cell(length(perms),1);
    S_next = cell(length(perms),1);
    H_next = zeros(length(perms),1);
    hL_next = zeros(length(perms),1);
    gap_next = zeros(length(perms),1);
    y_curr = cell(length(perms),1);
    parfor perm = perms
        
        % choose tie breaking order according to greatest gains of G or H w.r.t S
        if perm == 2 
            ties = G_marg_S;
        elseif perm == 3 
            ties = H_marg_S;
        else
            ties = randperm(n);
        end
        
        y = rho * v + greedy_algo_submodular(v, H.G, ties);
        
        % compute "prox" of f_L: argmin_{z in [0, 1]^n} f_L(z) + rho/2 ||z||^2 - z^Ty
%         if rho == 0
%             [H_approx, param_H_approx] = diff_sub_fct(F, param_F, @(S, param) sum(y(S)), struct);
%             S_next{perm} = submin_algo(H_approx, param_H_approx);
%             H_next(perm) = F(S_next{perm}, param_F) - G(S_next{perm}, param_G);
%             x_next{perm} = zeros(n,1);
%             x_next{perm}(S_next{perm}) = 1;
%         else 
            if warm_start % initialize sol of prox algo to current x
               x_init = x;
            else
               x_init = zeros(n,1);
            end
            [x_next{perm}, avg_obj_values, avg_gaps] = prox(H.F, V, y, rho, x_init);
            %fprintf("Average gap achieved at iter %d with perm %d: %f \n", iter, perm, avg_gaps(end))
            
%             plot(best_obj_values)
%             hold on
%             plot(avg_obj_values)
%             shg
%             if iter>1 % test if restart helps
%                 [x_next_norestart{perm}, obj_values_norestart] = prox(F, param_F, y, rho, zeros(n,1));
%                 figure
%                 plot(obj_values); hold on
%                 plot(obj_values_norestart);
%                 legend("prox-restart", "prox-norestart");
%                 keyboard
%             end
%             % test if solution returned by unrestricted prox of f_L projected
%             % on [0,1]^n is equivalent to restricted prox of f_L on
%             % [0,1]^n (not exactly but mnp sol has even better obj)
%             w = prox_operator_submodular(y/rho, 1/rho, F, param_F, 1e-5);
%             x_mnp = min(max(w,0),1);
%             mnp_prox_obj = 0.5*rho*norm(x_mnp)^2 - x_mnp'*y + lovasz_extension(x_mnp,F,param_F);
%             fprintf("mnp_prox_obj - pgm_prox_obj = %f \n", mnp_prox_obj - obj_values(end))

            % compute hL(x_next{perm}) and round 
            [w,values,order] = greedy_algo_submodular(x_next{perm},H.H);
            hL_next(perm) = w'*x_next{perm};
            [S_next{perm}, H_next(perm)] = round_lovasz(x_next{perm}, H.H, values, order);
            gap_next(perm) = avg_gaps(end);
            y_curr{perm} = y; % needed for criticality check
     end

    [H_next, perm_choices(iter)] = min(H_next);
    H_values(iter) = H_next;
    hL_values(iter) = hL_next(perm_choices(iter));
    x_next = x_next{perm_choices(iter)};
    S_next = S_next{perm_choices(iter)};
    gaps(iter) = gap_next(perm_choices(iter));
    iter_dist = norm(x - x_next); 
%     disp("Criticality of x_next:")
%     Tf_x_next = lovasz_extension(x, H.F) + 0.5*rho*norm(x)^2 - lovasz_extension(x_next, H.F) - 0.5*rho*norm(x_next)^2 - y_curr{ind}'*(x - x_next) % should be >=  - tol in prox and at convergence <= tol + tol in prox
%     Tg_x_next = lovasz_extension(x, H.G) + 0.5*rho*norm(x)^2 - lovasz_extension(x_next, H.G) - 0.5*rho*norm(x_next)^2 - y_curr{ind}'*(x - x_next) % should be <= 0 and at convergence >= - tol - tol in prox
%     if Tf_x_next < -avg_gaps(end)
%         keyboard
%     end
    if round
        x_next = zeros(n, 1);
        x_next(S_next) = 1;
        if iter>1
            obj_change = H_values(iter-1) - H_values(iter);
        else
            obj_change = - H_values(iter);
        end
%         disp("Criticality of rounded x_next:")
%         Tf_x_next = lovasz_extension(x, H.F) + 0.5*rho*norm(x)^2 - lovasz_extension(x_next, H.F) - 0.5*rho*norm(x_next)^2 - y_curr{ind}'*(x - x_next) % no bounds
%         Tg_x_next = lovasz_extension(x, H.G) + 0.5*rho*norm(x)^2 - lovasz_extension(x_next, H.G) - 0.5*rho*norm(x_next)^2 - y_curr{ind}'*(x - x_next) % should be <= 0
    else
        if iter>1
            obj_change = hL_values(iter-1) - hL_values(iter);
        else
            obj_change =  - hL_values(iter);
        end
    end
    
    x_prev = x;
    S = S_next;
    x = x_next;
    
    itertime(iter) = toc;
    if obj_change <= tol % iter_dist <= tol 
        if localmin
            [local_minimality, U] = check_localmin(H, S, V);
            if local_minimality>0
                fprintf("RegDC Algorithm converged after %d iterations but not to a local min, local minimality = %e \n", iter, local_minimality)
                S = U;
                x = zeros(n, 1);
                x(S) = 1;
                H_values(iter) = H_values(iter) - local_minimality;
                hL_values(iter) = H_values(iter);
            else
                fprintf("RegDC Algorithm converged after %d iterations to a local min, local minimality = %e \n", iter, local_minimality)
                break
            end
        else
            fprintf("RegDC Algorithm converged after %d iterations \n", iter)
            break
        end
    end
    
end
H_values = H_values(1:iter);
hL_values = hL_values(1:iter);
gaps = gaps(1:iter);
itertime = itertime(1:iter);
perm_choices = perm_choices(1:iter);

end