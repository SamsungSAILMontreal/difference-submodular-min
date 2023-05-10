classdef Entropy < SetFct
% Computes Entropy of a BINARY data set
    properties
        data = [];
        n_rows = 0;
    end
    methods
        function F = Entropy(data)
            F.data = uint64(data);
            F.n_rows = size(data,1);
        end
        
        function [val, F] = obj(F, A)
            %Returns the entropy wrt the columns (features) in A
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                if size(A,1)>1
                    A = A.';
                end
                if ~isempty(A)
                    A = sort(A);
                    counts = data_count(F.data(:, A), F.n_rows);
                    probs = counts / F.n_rows;
                    val = -1 * sum(probs .* log2(probs));

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
            % Returns the entropy wrt the columns (features) in A plus new
            % feature e
            if size(A,1)>1
                A = A.';
            end

            [val, F] = F.obj(A);
           if ismember(e, A)
              new_val = val;
           else
               cols = union(A,e);
               counts = data_count(F.data(:, cols), F.n_rows);
               probs = counts ./ F.n_rows;
               new_val = -1 * sum(probs .* log2(probs));

               F.current_set = cols;
               F.current_val = new_val; 
           end
           
        end
        
        function [new_val, F] = rmv(F, A, e)
            % Returns the entropy wrt the columns (features) in A without
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
                   counts = data_count(F.data(:, cols), F.n_rows);
                   probs = counts ./ F.n_rows;
                   new_val = -1 * sum(probs .* log2(probs));
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
