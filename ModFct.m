classdef ModFct < SetFct
% modular function F(A)= weights(A) where
% weights: non-negative weights of elements in V
    properties
        weights = 1;
    end
    methods
        function F = ModFct(weights)
            F.weights = weights;
        end
        
        function [val, F] = obj(F, A)
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                val = sum(F.weights(A));
                F.current_set = A;
                F.current_val = val;
            end
        end
        
        function [new_val, F] = add(F, A, e)
           [val, F] = F.obj(A);
           if ismember(e, A)
              new_val = val;
           else
               new_val = val + F.weights(e);
               F.current_set = union(A,e);
               F.current_val = new_val;
           end
        end
        
        function [new_val, F] = rmv(F, A, e)
           [val, F] = F.obj(A);
           if ~ismember(e, A)
              new_val = val;
           else
               new_val = val - F.weights(e);
               F.current_set = setdiff(A,e);
               F.current_val = new_val;
           end
        end
    end
end

