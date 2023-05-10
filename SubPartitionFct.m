classdef SubPartitionFct < SetFct
% submodular partition function F(A)= sum sqrt{c(A cap V_i} where
% V_1, ..., V_m is a partition of ground set V
% grp_edges: start and end indices of groups, with V_i = grp_edges(i):grp_edges(i+1)-1
% c: non-negative weights of elements in V
    properties
        c = 1;
        grp_edges = 1:2;
        m = 1;
        cA = 0;
    end
    methods
        function F = SubPartitionFct(c, grp_edges)
            F.c = c;
            F.grp_edges = grp_edges;
            F.m = length(grp_edges)-1;
        end
        
        function [val, F] = obj(F, A)
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                %%%% this is what changes in your example
                F.cA = zeros(F.m,1);
                for i=1:F.m
                   F.cA(i) = sum(F.c(intersect(A, F.grp_edges(i):F.grp_edges(i+1)-1)));
                end
                val = sum(sqrt(F.cA)); 
                %%%%%
                F.current_set = A;
                F.current_val = val;
            end
        end
        
        function [new_val, F] = add(F, A, e)
           [val, F] = F.obj(A);
           if ismember(e, A)
              new_val = val;
           else
               %%%% this is what changes in your example
               %find group Vi containing e
               ind = find(F.grp_edges <= e , 1, 'last');
               new_val = val + sqrt(F.cA(ind)+ F.c(e)) - sqrt(F.cA(ind));
               F.cA(ind) = F.cA(ind) + F.c(e);
               %%%%
               F.current_set = union(A,e);
               F.current_val = new_val; 
           end
        end
        
        function [new_val, F] = rmv(F, A, e)
           [val, F] = F.obj(A);
           if ~ismember(e, A)
              new_val = val;
           else
               %%%% this is what changes in your example
               %find group Vi containing e
               ind = find(F.grp_edges <= e , 1, 'last');
               new_val = val + sqrt(max(F.cA(ind)-F.c(e), 0)) - sqrt(F.cA(ind)); % added max in sqrt to avoid numerical errors  
               F.cA(ind) = F.cA(ind) - F.c(e);
               %%%% this is what changes in your example
               F.current_set = setdiff(A,e);
               F.current_val = new_val;
           end
        end
    end
end

