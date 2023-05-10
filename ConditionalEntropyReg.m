classdef ConditionalEntropyReg < SetFct
% Computes Conditional Entropy of a BINARY data set (data)
% with class c and a cardinality penalty weighted by lambda
    properties
        data = [];
        c = [];
        n_rows = 0;
        classEntropy = 0;
        lambda = 0;
    end
    methods
        function F = ConditionalEntropyReg(data, c, lambda)
            F.data = uint64(data);
            F.c = uint64(c);
            F.n_rows = size(data,1);
            F.lambda = lambda;
            F.classEntropy = - ((sum(c == 0)/length(c)) * log2((sum(c == 0)/length(c))) + (sum(c == 1)/length(c)) * log2((sum(c == 1)/length(c))));
        end
        
        function [val, F] = obj(F, A)
            %Returns the conditional entropy wrt the columns (features) in A
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                if size(A,1)>1
                    A = A.';
                end
                if ~isempty(A)
                    A = sort(A);
                    counts = data_count([F.data(:,A),F.c],F.n_rows);
                    probs = counts ./ F.n_rows;
                    ent = -1 * sum(probs .* log2(probs));
                    val = ent - F.classEntropy + F.lambda * length(A);

                    F.current_set = A;
                    F.current_val = val;  
                else
                    val = 0;
                    F.current_set = A;
                    F.current_val = val;
                end
            end
        end
        
        function [new_val, F] = add(F, A, e)
            % Returns the conditional entropy wrt the columns (features) in A plus new
            % feature e
            if size(A,1)>1
                A = A.';
            end

            [val, F] = F.obj(A);
           if ismember(e, A)
              new_val = val;
           else
               cols = union(A,e);
               counts = data_count([F.data(:,cols),F.c],F.n_rows);
               probs = counts ./ F.n_rows;
               ent = -1 * sum(probs .* log2(probs));
               new_val = ent - F.classEntropy + F.lambda * length(cols);
  
               F.current_set = cols;
               F.current_val = new_val; 
           end
        end
        
        
        function [new_val, F] = rmv(F, A, e)
            % Returns the conditional entropy wrt the columns (features) in A without
            % feature e
            if size(A,1)>1
                A = A.';
            end
            
            [val, F] = F.obj(A);
            if ~ismember(e, A)
              new_val = val;
            else
               cols = setdiff(A,e);
               if ~isempty(cols)
                   counts = data_count([F.data(:,cols),F.c],F.n_rows);
                   probs = counts ./ F.n_rows;
                   ent = -1 * sum(probs .* log2(probs));
                   new_val = ent - F.classEntropy + F.lambda * length(cols);
                   F.current_set = cols;
                   F.current_val = new_val;
               else
                   new_val = 0;
                   F.current_set = cols;
                   F.current_val = new_val;
               end
           end
        end
   end
end

