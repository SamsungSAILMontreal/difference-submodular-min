function [S, x, H_values, hL_values, gaps, prox_max_gaps, perm_choices, itertime, FW_niters] = regCompleteDC(H, V, rho, perms, prox, maxiter, tol, FW_maxiter, FW_gap, warm_start, ls, round, localmin)
% Compute  min_{S} H(S):= F(S) - G(S)
% by applying complete DC algo to the continuous problem
% min_x f_L(x) + delta_[0,1]^n(x) + rho ||x||^2 - g_L(x) - rho ||x||^2
% where f_L and g_L are the Lovasz extensions of F and G
% INPUT
% H: difference of submodular functions H.H(S) = H.F(S)- H.G(S), with
% H.H, H.F and H.G of class SetFct
% V: ground set
% rho: regularization parameter
% prox: proximal operator of f_L restricted to [0,1]^n returns
% [w, obj_values] = prox(F,param_F,y,rho, w_init) where 
% w = argmin_{w in [0, 1]^n} f_L(w) + rho/2 ||w||^2 - w^Ty
% maxiter: maximum number of iterations
% localmin: use local minimality as additional stopping criterion

% OUTPUT
% S: argmin H(S)
% Hvalues: H(S) at every iteration

n = length(V);
x = zeros(n,1);
S = [];
% y should be in rho x + subgrad g_L(x)
y = rho * x + greedy_algo_submodular(x, H.G);
H_values = zeros(maxiter, 1);
hL_values = zeros(maxiter, 1);
gaps = zeros(maxiter, 1);
prox_max_gaps = zeros(maxiter, 1);
perm_choices = zeros(maxiter, 1);
itertime = zeros(maxiter, 1);
FW_niters = zeros(maxiter, 1);
H_marg_S = zeros(n,1);
G_marg_S = zeros(n,1);

if rho==0 %can't use line search if rho=0
    ls = false;
end

function [obj, grad, info] = h_dual_approx(y, x, info)
    % compute <y,x> - f*(y) and its gradient where
    % f* = (f_L + rho/2 ||.||^2 + delta_[0,1]^n)*
    % f* = max_{z in [0, 1]^n} z^Ty - f_L(z) - rho/2 ||z||^2 
    % grad f* = argmin_{z in [0, 1]^n} f_L(z) + rho/2 ||z||^2 - z^Ty
    
    if warm_start
        z_init = info.prox_f_init;
    else
        z_init = zeros(n,1);
    end
    [info.prox_f_init, prox_values, prox_gaps] = prox(H.F,V, y, rho, z_init);
    obj = y'*x + prox_values(end);
    grad = x - info.prox_f_init;
    info.prox_gap = max(prox_gaps(end), info.prox_gap);
end

function v = LO(w, x)
    % compute argmin_v {w^Tv : v \in rho x + subgrad g_L(x))}
    %         = argmax_s {-w^Ts: s \in subgrad g_L(x)} + rho x
    % sort indices in decreasing order of entries in x, break ties
    % according to decreasing order of entries in w
    v = rho * x + greedy_algo_submodular(x, H.G, -w);
end

delete(gcp('nocreate')); % ensure no parallel pool is running
ncpus = str2num(getenv('SLURM_CPUS_ON_NODE'))
if isempty(ncpus) % running on laptop
    parpool(3);
else
    parpool(min(str2num(getenv('SLURM_CPUS_ON_NODE')), 3));
end
% Start the clock.
tic;
iter = 0;
prox_nops = 0;
while prox_nops <  maxiter
    iter = iter + 1;
    % y in argmin_y {<y, x> - f*(y): y in rho x + subgrad g_L(x)} where 
    % f* = (f_L + rho/2 ||.||^2 + delta_[0,1]^n)*
    h_dual_approx_x = @(y, prox_f_prev) h_dual_approx(y, x, prox_f_prev);
    if warm_start
        info.prox_f_init = x;
    end
    info.prox_gap = 0;
    
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
    
    % it's not possible to warm start FW, since the constraint
    % rho x + subgrad g_L(x) changes at each iteration. Instead we init
    % with best feasible point among ones obtained by each perm
    y_init = cell(length(perms),1);
    prox_f = cell(length(perms),1);
    FW_obj_init = zeros(length(perms),1);
    parfor perm = perms
        % choose tie breaking order according to greatest gains of G or H w.r.t S
        if perm == 2 
            ties = G_marg_S;
        elseif perm == 3 
            ties = H_marg_S;
        else
            ties = randperm(n);
        end        
        y_init{perm} = rho * x + greedy_algo_submodular(x, H.G, ties);
        [FW_obj_init(perm), ~, info_perm] = h_dual_approx_x(y_init{perm}, info)
        prox_f{perm} = info_perm.prox_f_init
    end
    [FW_obj_init, perm_choices(iter)] = min(FW_obj_init);
    y_init = y_init{perm_choices(iter)};
    if warm_start
        info.prox_f_init = prox_f{perm_choices(iter)};
    end    
    LO_x = @(w) LO(w, x);
    
    % Lipschitz constant of gradient of <y, x> - f*(y) is 1/rho
    [y, FW_obj_values, FW_gaps, info] = FW(h_dual_approx_x, LO_x, 1/rho, min(FW_maxiter, maxiter - prox_nops), FW_gap, false, true, y_init, info);
    
    if warm_start % initialize sol of prox algo with last prox computed in h_dual_approx
       x_init = info.prox_f_init;
    else
       x_init = zeros(n,1);
    end
    [x_next, prox_values, prox_gaps] = prox(H.F, V, y, rho, x_init);
    info.prox_gap = max(prox_gaps(end), info.prox_gap);
    
    % compute hL(x_next) and round 
    [w,values,order] = greedy_algo_submodular(x_next,H.H);
    hL_values(iter) = w'*x_next;
    [S, H_values(iter)] = round_lovasz(x_next, H.H, values, order);
    gaps(iter) = FW_gaps(end);
    prox_max_gaps(iter) = info.prox_gap; % max gap of all prox operations done this iter
    FW_niters(iter) = length(FW_gaps);
    iter_dist = norm(x - x_next); 
    if round
        x_next = zeros(n, 1);
        x_next(S) = 1;
        if iter>1
            obj_change = H_values(iter-1) - H_values(iter);
        else
            obj_change = - H_values(iter);
        end
    else
        if iter>1
            obj_change = hL_values(iter-1) - hL_values(iter);
        else
            obj_change =  - hL_values(iter);
        end
    end
    
    x = x_next;
    prox_nops = sum(FW_niters) + iter;
    itertime(iter) = toc;
    if obj_change <= tol %iter_dist <= tol 
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
prox_max_gaps = prox_max_gaps(1:iter);
itertime = itertime(1:iter);
perm_choices = perm_choices(1:iter);
FW_niters = FW_niters(1:iter);
end