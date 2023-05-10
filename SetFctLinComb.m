classdef SetFctLinComb < SetFct
    % Creates a linear combination of set functions F = sum_i Fs{i}*weights(i)
    
    properties
        Fs = {};
        weights = 1;
    end
    
    methods
        function F = SetFctLinComb(Fs, weights)
            F.Fs = Fs;
            F.weights = weights;
        end
        
        function [val, F] = obj(F, A)
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                val = 0;
                for i=1:length(F.Fs)
                    [vali, F.Fs{i}] = obj(F.Fs{i}, A);
                    val = val + F.weights(i)*vali;
                end
                F.current_set = A;
                F.current_val = val;
            end
        end
        
        function [new_val, F] = add(F, A, e)
           [val, F] = F.obj(A);
           if ismember(e, A)
              new_val = val;
           else
                new_val = 0;
                for i=1:length(F.Fs)
                    [new_vali, F.Fs{i}] = add(F.Fs{i}, A, e);
                    new_val = new_val + F.weights(i)*new_vali;
                end
                F.current_set = union(A,e);
                F.current_val = new_val;
           end
        end
        
       function [new_val, F] = rmv(F, A, e)
           [val, F] = F.obj(A);
           if ~ismember(e, A)
              new_val = val;
           else
                new_val = 0;
                for i=1:length(F.Fs)
                    [new_vali, F.Fs{i}] = rmv(F.Fs{i}, A, e);
                    new_val = new_val + F.weights(i)*new_vali;
                end
                F.current_set = setdiff(A,e);
                F.current_val = new_val;
           end
        end
    end
end

