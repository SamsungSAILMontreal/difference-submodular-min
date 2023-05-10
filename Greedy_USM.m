function [X_max, fvalues, itertime] = Greedy_USM(f, V, randomized)
% Unconstrained submodular maximization: max_S F(S)
% f: submodular function of class SetFct
% V: ground set
% if randomized = 0
% We use the "Deterministic USM" algo from SMunconstrained(journal)
% which achieves a 1/3 approximation
% if randomized = 1
% We use the "Randomized USM" algo from SMunconstrained(journal)
% which achieves a 1/2 approximation in expectation
% if randomized = 2 (TODO: modify code to use SetFct interface for this option)
% We use the "Deterministic Unconstrained" algorithm from 
% "Deterministic Algorithms for Submodular Maximization Problems"
% which achieves 1/2 approximation
n = length(V);
fvalues = zeros(n,1);
itertime = zeros(n,1);

if nargout>=3
    % Start the clock.
   tic;
end

if randomized ~= 2
    X = [];
    Y = V;
    [fX_val, fX_obj] = f(X);
    [fY_val, fY_obj] = f(Y);
    
       for i = V    
            [fXi_val, fXi_obj] = add(fX_obj, X, i); 
            f_marg_X = fXi_val - fX_val; %f_marg(X, i);
            [fYmi_val, fYmi_obj] = rmv(fY_obj, Y, i); 
            f_marg_Y = fYmi_val - fY_val; %-f_marg(setdiff(Y,i),i);

            if randomized
                a = max(f_marg_X,0);
                b = max(f_marg_Y,0);
                if a == 0 && b == 0
                    p = 1;
                else
                    p = a/(a + b);
                end

                if binornd(1,p)
                    X = union(X, i);
                    fX_val = fXi_val;
                    fX_obj = fXi_obj;
                else
                    Y = setdiff(Y,i);
                    fY_val = fYmi_val;
                    fY_obj = fYmi_obj;
                end
            else

                if  f_marg_X >= f_marg_Y
                    X = union(X, i);
                    fX_val = fXi_val;
                    fX_obj = fXi_obj;
                else
                    Y = setdiff(Y,i);
                    fY_val = fYmi_val;
                    fY_obj = fYmi_obj;
                end
            end
            fvalues(i) = max(fX_val, fY_val);
            if nargout>=3
                itertime(i) = toc;
            end
       end
       
       X_max = X;
       
else
    m = 1;
    p = 1;
    X = {[]};
    Y = {V};
    
    for i = 1:n
        
        a = zeros(m,1);
        b = zeros(m,1);

        for j = 1:m
            a(j) = f_marg(X{j},i);
            b(j) = -f_marg(setdiff(Y{j},i),i);
        end
        
        %% find extreme point of (P) by solving fractional knapsack problem (see Claim A.1)
        v = p.*(a - 3*b);
        s = p.*(b - 3*a);
        B = p'*(b - 2*a);
        
%         cvx_begin 
%             cvx_precision best
%             variable z_cvx(m);
%             maximize (z_cvx'*v)
%             subject to 
%                 z_cvx'*s <= B
%                 z_cvx <= 1
%                 z_cvx >= 0
%         cvx_end
        
        % take all items with vj ≥ 0 and sj ≤ 0 & omit all items of vj < 0 and sj ≥ 0
        z = zeros(m,1);
        z((v >= 0) & (s <= 0))=1;
        % update budget
        B = B - z'*s;
        
        vs = v./s;
        items = 1:m;
        
        % sort the positive items (vj ≥ 0 and sj > 0) in non-increasing order of vj/sj
        pos_items = items((v >= 0) & (s > 0));
        [~,sorted_pos] = sort(vs(pos_items),'descend');
        
        % sort the negative items (vj ≤ 0 and sj < 0) in non-decreasing order of vj/sj
        neg_items =  items((v <= 0) & (s < 0));
        [~,sorted_neg] = sort(vs(neg_items),'ascend');
        
        j_pos = 1;
        while B> 0 && j_pos<= numel(pos_items)
         % add positive items until B reaches 0 (or we are out of positive items)
            ind = pos_items(sorted_pos(j_pos)); 
            z(ind) = min(B/s(ind),1);
            if z(ind) == 1
               j_pos = j_pos + 1;
            end
            B = B - z(ind)*s(ind);
        end
        
        j_neg = 1;
        while B<0
         % add negative items until B reaches 0 
            ind = neg_items(sorted_neg(j_neg)); 
            z(ind) = min(B/s(ind),1);
            if z(ind) == 1
               j_neg = j_neg + 1;
            end
            B = B - z(ind)*s(ind);
        end
        
        while vs(j_pos)>= vs(j_neg) && j_pos<= numel(pos_items) && j_neg<= numel(neg_items)
          % add positive & negative items in such a way that B stays 0
            ind_pos = pos_items(sorted_pos(j_pos)); 
            ind_neg = neg_items(sorted_neg(j_neg)); 
            
            z(ind_pos) = min(z(ind_pos) - s(ind_neg)/s(ind_pos),1);
            if z(ind_pos) == 1
                j_pos = j_pos + 1;
            end
                
            z(ind_neg) =  z(ind_neg) - z(ind_pos)*s(ind_pos)/s(ind_neg);
            if z(ind_neg) == 1
                j_neg = j_neg + 1;
            end
        end
        
        w = 1 - z;
        %% update distribution
        m_old = m;
        for j = 1:m_old
            if z(j) == 1
                X{j} = [X{j},i];
            elseif w(j) ==1
                Y{j} = setdiff(Y{j},i);
            else %fractional sol 
                p(m+1,1) = w(j)*p(j);
                X{m+1} = X{j};
                Y{m+1} = setdiff(Y{j},i);
                
                p(j) = z(j)*p(j);
                X{j} = [X{j},i];
                
                m = m+1;
            end
        end 
    end
    
    %% return set in the support of the distribution with largest f value
    f_max = 0;
    X_max = [];
    
    for j=1:m
        if f(X{j})> f_max
            f_max = f(X{j});
            X_max = X{j};
        end
    end
    
end