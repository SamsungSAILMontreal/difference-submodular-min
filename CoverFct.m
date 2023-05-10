classdef CoverFct < SetFct
% Cover function F(A)= | U_{i in A} G_i|^alpha where
% W is a sparse matrix with W(i, j) = 1 if j in G_i, 0 otherwise
% alpha is any power <= 1

    properties
        W = sparse(1,1);
        alpha = 1;
        coverA = sparse(1,1);
    end
    methods
        function F = CoverFct(W, alpha)
            F.W = W;
            F.alpha = alpha;
        end
        
        function [val, F] = obj(F, A)
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                F.coverA = any(F.W(A,:), 1);
                val = full(sum(F.coverA))^F.alpha;
                F.current_set = A;
                F.current_val = val;
            end
        end
        
        function [new_val, F] = add(F, A, e)
           [val, F] = F.obj(A);
           if ismember(e, A)
              new_val = val;
           else
               F.coverA = F.W(e,:)~=0 | F.coverA;
               new_val = full(sum(F.coverA))^F.alpha;
               F.current_set = union(A,e);
               F.current_val = new_val;
           end
        end
        % use rmv from base class, there's no more efficient implementation 
    end
end

