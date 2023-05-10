function [x, obj_values, gap_values, info] = FW(f, LO, L, maxiter, gap_tol, ls, concave, x_init, info)
%% Solve min_{x in C} f(x) using Frank-Wolfe algorithm
% Input:
% f: smooth, L-Lip gradient fct, f(x, info) = [f value at x, grad of f at x, info]
% LO: linear minimization oracle, LO(z) = armin_{x in C} <z, x>
% for step size choice see for example 
% http://fa.bianp.net/blog/2018/notes-on-the-frank-wolfe-algorithm-part-i/
% for concave functions, we can use the greedy step size gamma = 1, and the
% obj is almost monotonically decreasing (up to error of grad)
obj_values = zeros(1, maxiter);
gap_values = zeros(1, maxiter);
x = x_init;
for iter = 1:maxiter
    x_prev = x;
    [fx, gradfx, info] = f(x, info);
    v = LO(gradfx);
    d = v - x;
    obj_values(iter) = fx;
    gap_values(iter) = - gradfx'*d;
%     if  gap_values(iter)<0
%         keyboard
%     end
    if gap_values(iter) < gap_tol
        fprintf('FW reached small duality gap (gap=%g) after %d iterations \n', gap_values(iter), iter);       
        break;
    end
    if concave
        gamma = 1;
    else
        if ls
            gamma = min(gap_values(iter)/(L*norm(d)^2), 1);
        else
            gamma = 2/(iter+1);
        end
    end
    x = (1 - gamma)*x_prev + gamma*v;
end
obj_values = obj_values(1:iter);
gap_values = gap_values(1:iter);
end