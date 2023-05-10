function [x_primal,x_dual,dual_values,primal_values,gaps,added1,added2,itertime] = minimize_submodular_FW_minnormpoint(F,V,maxiter,sfm,gap,optimal_primal,starting_point)
% Compute  min_x .5 * || x ||^2 + f(x) using Frank-Wolfe on the dual (using
% the min norm point algorithm)
%
% if sfm=1, perform submodular function minimization (and adapt gaps and
% function values accordingly)
%
%
% INPUT
% F: submodular function of class SetFct
% maxiter: maximum number of iterations
% sfm: 1 if minimizing submodular function, 0 otherwise
%
% OUTPUT
% x: argmin. If sfm=1, then outputs the subset
% values of cost function
% gaps: certified optimality gaps at each iteration
if nargin<6, optimal_primal = -Inf; end
if nargin<4, sfm=0; end

if nargout>=8
    tic;
end

n = length(V);

display=1;

 
if nargin < 7
    % use random initialization
    direction = rand(n,1);
     
    x = greedy_algo_submodular(direction,F);
else
    x = starting_point;
end
w = 1;%M: these are the coefficients associated with the chosen atoms
X = x;%M: we initialize with an atom
nstep=1;
iter = 0;
bestvalue = Inf;


while iter < maxiter
    iter = iter + 1;
    
    % step 1
    % a
    x = X*w;

    direction = - X*w;

    % b %M: we call the LO with - gradient,i.e., -x
    [xx,Fvalues,order] = greedy_algo_submodular(direction,F);

    % c


    if sfm
        % output what's needed for SFM
        [a,b] = min(Fvalues);%M: Fvalues are the sublevels of x, {x<=theta}s
        if min(a,0) < bestvalue
            if a < 0
                % allow empty set to be optimal
                Aopt = order(1:b);
                bestvalue = a;
            else
                bestvalue = 0;
                Aopt = [];
            end
        end
        
        if nargout>=8
            itertime(iter) = toc;
        end
        dual_values(iter) = sum(min(x,0));
        primal_values(iter) = min(a,0);
        
        % compute traditional value
        [added1(iter), F] = F(find(x<-1e-15));
        [added2(iter), F] = F(find(x<+1e-15));
       
        gaps(iter) = a - dual_values(iter);

        if gaps(iter) < gap
            if display, fprintf('stop because reached duality gap on SFM %f\n',gaps(iter)); end
            break;
        end

        ww = pav( -xx(flipud(order)) );
        ww(flipud(order)) = ww;
        gap_quad = ww'*xx + .5 * ww' * ww + .5 * x' * x;%M: primal-dual gap
%         if gap_quad  < gap.^2 / n
%             % stop!
%             if display, fprintf('stopped at step 1c - reached upper bound on duality gap \n'); end
%             break;
%         end
%         
        if primal_values(iter) < optimal_primal
            if display, fprintf('stop because reached required accuracy on SFM %f\n',optimal_primal - primal_values(iter)); end
            break;
        end

    else



        % compute values of the function
        dual_values(iter) = - .5 * x' * x;
        ww = pav( -xx(flipud(order)) );
        ww(flipud(order)) = ww;
        
        if nargout>=8
            itertime(iter) = toc;
        end
        added1(iter)  = xx'* ( - x) + .5 * x'*x; % compute value without PAV
        primal_values(iter)  = ww'*xx + .5 * ww' * ww;
         
        gaps(iter) = primal_values(iter) - dual_values(iter);
        added2(iter)  = added1(iter) - dual_values(iter); % compute gap without PAV
        if gaps(iter) < gap
            % stop!
            if display, fprintf('stopped at step 1c - reached upper bound on duality gap \n'); end
            break;
        end

    end


    % d
    if min(sum( (X - repmat(xx,1,size(X,2))).^2 , 1 ) ) < 1e-12 *max(sum(X.^2,1));
        % stop!
        if display, fprintf('stopped at step 1d\n'); end
        break;
    end

    % e
    X = [X, xx];%M: add the new atom u got from LO to the atoms set (what if it's already there?)
    w = [w; 0];

    iterloc =1;
    while 1,
        iterloc = iterloc + 1;
        if iterloc > 100*n,
            fprintf('probably looping between step 2 and 3, exit \n');
            iterloc = 0;
            break;
        end
        % step 2
        % a
        try
            R = chol( X'*X + 1/size(X,2) * norm(X,'fro').^2 + 1e-15 * 1/size(X,2) * norm(X,'fro').^2 * eye( size(X,2)) );
        catch
            % not positive definite
            iterloc = 0;
            fprintf('not positive definite when adding new point in step 2, exit\n');
            break;
        end
        v = R \ ( R' \ ones(size(X,2),1) );
        v = v / sum(v);

        % checking
        % X'* ( X * v );

        % b
        if all( v > 1e-12 ),
            w = v;
            break;
            % go to step 1
        else

            % step 3
            % a
            ind = find( w - v > 1e-12 );%these are the atoms that still have weight >0

            if isempty(ind)%if there's none get out
                iterloc = 0;
                fprintf('can''t do line seaerch in step 2, exit\n');
                break;
            end
            %b
            theta = min( 1, min ( w(ind)./(w(ind)-v(ind))));
            %c
            w = theta * v + (1-theta) * w;
            %d-e
            torem = find( w < 1e-12);
            w(torem) = [];
            X(:,torem) = [];
            w = w / sum(w);
            % go to step 2
        end
    end
    if iterloc==0, break; end
end

if sfm
    x_dual = x;
    x_primal = Aopt; %find( -x_dual >= 1e-8 );
else
    x_primal = ww;
    x_dual = x;
end
% if nargout>=8
%     itertime(iter+1) = toc;
% end
% 
% [xx,Fvalues] = greedy_algo_submodular(-x,F,param_F);
% 
% % compute values of the function
% if sfm
% 
%     [a,b] = min(Fvalues);
%     if min(a,0) < bestvalue
%         if a < 0
%             % allow empty set to be optimal
%             Aopt = order(1:b);
%             bestvalue = a;
%         else
%             bestvalue = 0;
%             Aopt = [];
%         end
%     end
%     dual_values(iter+1) = sum(min(x,0));
%     primal_values(iter+1) = min(a,0);
%     gaps(iter+1) = primal_values(iter+1) - dual_values(iter+1);
%     x_dual = x;
%     x_primal = Aopt; %find( -x_dual >= 1e-8 );
% 
% else
% 
%     dual_values(iter+1) = - .5 * x' * x;
%     ww = pav( -xx(flipud(order)) );
%     ww(flipud(order)) = ww;
%     added1(iter+1)  = xx'* ( - x) + .5 * x'*x;
%     primal_values(iter+1)  = ww'*xx + .5 * ww' * ww;
%     gaps(iter+1) = primal_values(iter+1) - dual_values(iter+1);
%     added2(iter+1)  = added1(iter+1) - dual_values(iter+1);
%     x_primal = ww;
%     x_dual = x;
% 
% end

