function [w_hat, obj_values, avg_gaps, itertime, w_hat_round, obj_values_round] = prox_lovasz_hypercube(F,V,y,rho, maxiter, tol, step_size, w_init, L, round)
% Compute "prox" of f_L: argmin_{w in [0, 1]^n} f_L(w) + rho/2 ||w||^2 - w^Ty
% using projected subgradient descent
% See section 3.2.3 p. 132 in "Introductory Lectures on Convex Programming"
% and section 3.1 p. 24 in "Theory of Convex Optimization for Machine
% Learning" for non-strongly convex case (rho=0)
% See "A simpler approach to obtaining an O(1/t) convergence rate for the 
% projected stochastic subgradient method" and section 3.4.1 p. 37 in 
% "Theory of Convex Optimization for Machine Learning" for strongly convex 
% case (rho > 0)
% F submodular function of class SetFct
% L Lipschitz constant of f_L set to F(V) if F is monotone, otherwise to 3*max_S |F(S)| if known, and Inf otherwise 
% round return rounded solution based on F 
n = length(V);
if nargin < 7
    w = zeros(n,1);
else 
    w = w_init; % assumed to be integral if round
end
w_avg = w;
%w_best = w;
best_obj = 0.5*rho*norm(w)^2 - w'*y + lovasz_extension(w,F);
best_obj_round = best_obj;
s_avg = zeros(n,1); % I tried warm starting s too, but it didn't help much 
D = sqrt(n);
if L == Inf % use upper bound sqrt(sum_i F(i)^2)
    L = 0;
    for i = V
        L = L + F(i)^2;
    end
    L = sqrt(L); 
end
L = L + norm(y)+ rho* D;

itertime = zeros(maxiter, 1);
obj_values_round = zeros(maxiter, 1);
obj_values = zeros(maxiter, 1);
avg_gaps = zeros(maxiter, 1); % avg gap is guaranteed to converge (see Bach_learning_new section 7.2). 
%I also computed gaps of current iterates and best iterates, but these did not converge in the few exps I ran 

% Start the clock.
   tic;
   
for iter =1:maxiter
    [s,values,order] = greedy_algo_submodular(w,F);
    s = s - y;
    
    subgrad = rho * w + s;
    
    if rho==0 
        w_avg = (w_avg * (iter-1) + w)/iter;
        s_avg = (s_avg * (iter-1) + s)/iter;
    else
        w_avg = (w_avg * (iter-1) + 2*w)/(iter+1);
        s_avg = (s_avg * (iter-1) + 2*s)/(iter+1);
    end
    
    if round % return sol with best obj (to be consistent with what we did in "Optimal approximation for unconstrained non-submodular minimization")
        [w_round, F_round] = round_lovasz(w, F, values, order);
        obj_values_round(iter) = 0.5*rho*length(w_round)^2 - sum(y(w_round)) + F_round;
        if obj_values_round(iter) < best_obj_round
            w_hat = w;
            w_hat_round = w_round;
            obj_values(iter) = 0.5*rho*norm(w)^2 - w'*y + s'*w;
            best_obj = obj_values(iter);
            best_obj_round = obj_values_round(iter);
        else
            obj_values_round(iter) = best_obj_round;
            obj_values(iter) = best_obj;
        end
    else 
        w_hat = w_avg;
        obj_values(iter) = 0.5*rho*norm(w_avg)^2 - w_avg'*y + lovasz_extension(w_avg,F);
    end

    
    % compute duality gap
    z_avg = min(max(-s_avg/rho, 0),1);
    avg_gaps(iter) = obj_values(iter) - s_avg'*z_avg - 0.5*rho*norm(z_avg)^2;
   
    itertime(iter) = toc;
    if avg_gaps(iter)<= tol  %(iter > 1 && abs(avg_obj_values(iter - 1) - avg_obj_values(iter)) <= tol) 
        fprintf("prox algorithm converged after %d iterations \n", iter);       
        break;
    end

    if rho ==0 
        eta = D/(L*sqrt(iter));
    else 
        % In "A simpler approach to obtaining an O(1/t) convergence rate
        % for the projected stochastic subgradient method" they found that
        % using 1/rho*t step size is better than 2/rho*(iter+1)
        if step_size==1
            eta = 1/(rho*iter);
        else
            eta = 2/(rho*(iter+1));
        end
    end
    w = w - eta * subgrad;
    w = min(max(w,0),1);

end
obj_values_round = obj_values_round(1:iter);
obj_values = obj_values(1:iter);
avg_gaps = avg_gaps(1:iter);
itertime = itertime(1:iter);
end