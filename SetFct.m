classdef (Abstract) SetFct 
% Base class for set functions
% keeps track of current state of function to allow for efficient
% implementation of F(A U e) and F(A \ e) after evaluting F(A)  
% This class is based on sfo_fn class in SFO Matlab Toolbox by Andreas Krause (krausea@gmail.com), 
% available at: https://www.mathworks.com/matlabcentral/fileexchange/20504-submodular-function-optimization
% but with slightly different functionality (mainly add, rmv functions
% return updated function, not just value) and with new class definition
% syntax of Matlab

   properties
      current_set = -1;
      current_val = 0;
   end
   methods (Abstract)
       % obj should evaluate F(A) and update current_set and current_val
       [val, F] = obj(F, A)
   end
   methods
       function [new_val, F] = add(F, A, e)
           [new_val, F] = F.obj(union(A,e));
       end
       function [new_val, F] = rmv(F, A, e) 
           [new_val, F] = F.obj(setdiff(A,e));
       end
       function [val, F] = subsref(F,S)
           switch S(1).type
            case '()'
              A = S(1).subs{:};
              [val, F] = F.obj(A);
            case '{}'
               error('{} indexing not supported');
            case '.'  
               % this allows accessing properties like F.current_set, 
               % but accessing methods as F.add will only return val, use add(F, A, e) instead
               val = builtin('subsref', F, S);  
           end
       end
   end
end


