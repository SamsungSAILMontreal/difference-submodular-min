classdef MutualInformation < SetFct
% Computes Mutual Information of a BINARY data set (data)
% with class c
    properties
        count_dict = containers.Map('KeyType','char','ValueType','double');
        index_dict = containers.Map('KeyType','char', 'ValueType','any');
        data = [];
        c = [];
        n_rows = 0;
        classEntropy = 0;
    end
    methods
        function F = MutualInformation(data, c)
            F.data = data;
            F.c = c;
            F.n_rows = size(data,1);
            F.classEntropy = - ((sum(c == 0)/length(c)) * log2((sum(c == 0)/length(c))) + (sum(c == 1)/length(c)) * log2((sum(c == 1)/length(c))));
        end
        
        function [val, F] = obj(F, A)
            if size(A,1)>1
                A = A.';
            end
            new_count_dict = containers.Map('KeyType','char','ValueType','double');
            new_index_dict = containers.Map('KeyType','char', 'ValueType','any');
            
            new_cols = sort(A);
            if isequal(A, F.current_set)
                val =  F.current_val;
            else
                for i = 1:F.n_rows
                    perm = F.data(i,new_cols);
                    key_ent = num2str(perm);
                    key_joint = num2str([perm,F.c(i)]);
                    %Try to update mapping using key, if key is not present in map, it is a
                    %newly observed permutation, so map(key) = 1.
                    try
                        new_count_dict(key_ent) = new_count_dict(key_ent) + 1;
                        new_index_dict(key_ent) = [new_index_dict(key_ent), i];
                    catch
                        new_count_dict(key_ent) = 1;
                        new_index_dict(key_ent) = i;
                    end
                    
                    try
                        new_count_dict(key_joint) = new_count_dict(key_joint) + 1;
                        new_index_dict(key_joint) = [new_index_dict(key_joint), i];
                    catch
                        new_count_dict(key_joint) = 1;
                        new_index_dict(key_joint) = i;
                    end
                end
                
                F.count_dict = new_count_dict;
                F.index_dict = new_index_dict;
                
                ent = 0;
                joint = 0;
                key_set = keys(F.count_dict);
                m = length(A);

                for j = 1:length(key_set)
                    key = char(key_set(j));
                    key_array = str2num(key);
                    if length(key_array) == m
                        prob = F.count_dict(key)/F.n_rows;
                        ent = ent - (prob) * log2(prob);
                    else
                        prob = F.count_dict(key)/F.n_rows;
                        joint = joint - (prob) * log2(prob);
                    end
                end
                val = -1*(ent - joint) - F.classEntropy;
                F.current_set = new_cols;
                F.current_val = val;
            end
        end
        
        function [new_val, F] = add(F, A, e)
            if size(A,1)>1
                A = A.';
            end
            new_count_dict = containers.Map('KeyType','char','ValueType','double');
            new_index_dict = containers.Map('KeyType','char', 'ValueType','any');
            n = length(A);
            [val, F] = F.obj(A);
            if ismember(e, A)
                new_val = val;
            else
               key_set = keys(F.count_dict);
               cols = sort([F.current_set, e]);
                for i = 1:length(key_set)
                   key_string = char(key_set(i));
                   current_key = str2num(key_string);
                   if length(current_key) == n
                       current = F.current_set;
                   else
                       current = [F.current_set, Inf];
                   end
                   upper = current_key(current > e);
                   lower = current_key(current < e);
                   rows = F.index_dict(char(key_set(i)));
                   key_to_add_0 = num2str([lower, 0, upper]);
                   key_to_add_1 = num2str([lower, 1, upper]);
                    for j = 1:length(rows)
                        if F.data(rows(j), e) == 0
                            try
                                new_count_dict(key_to_add_0) = new_count_dict(key_to_add_0) + 1;
                                new_index_dict(key_to_add_0) = [new_index_dict(key_to_add_0), rows(j)];
                            catch
                                new_count_dict(key_to_add_0) = 1;
                                new_index_dict(key_to_add_0) = rows(j);
                            end
                        end
                    end
                    
                    try
                        temp = F.count_dict(key_string) - new_count_dict(key_to_add_0);
                        if temp ~= 0
                            new_count_dict(key_to_add_1) = F.count_dict(key_string) - new_count_dict(key_to_add_0);
                            new_index_dict(key_to_add_1) = setdiff(F.index_dict(key_string), new_index_dict(key_to_add_0));
                        end
                    catch
                        new_count_dict(key_to_add_1) = F.count_dict(key_string);
                        new_index_dict(key_to_add_1) = F.index_dict(key_string);
                    end
                end
                
                F.count_dict = new_count_dict;
                F.index_dict = new_index_dict;
               
                ent = 0;
                joint = 0;
                key_set = keys(F.count_dict);
                m = n + 1;

                for j = 1:length(key_set)
                    key = char(key_set(j));
                    key_array = str2num(key);
                    if length(key_array) == m
                        prob = F.count_dict(key)/F.n_rows;
                        if prob ~= 0
                            ent = ent - (prob) * log2(prob);
                        end
                    else
                        prob = F.count_dict(key)/F.n_rows;
                        if prob ~= 0
                            joint = joint - (prob) * log2(prob);
                        end
                    end
                end
                new_val = -1*(ent - joint) - F.classEntropy;
                
                F.current_set = cols;
                F.current_val = new_val;
               
           end
           
        end
        
        function [new_val, F] = rmv(F, A, e)
            if size(A,1)>1
                A = A.';
            end
            new_count_dict = containers.Map('KeyType','char','ValueType','double');
            new_index_dict = containers.Map('KeyType','char', 'ValueType','any');
            n = length(A);
            [val, F] = F.obj(A);
            if ~ismember(e, A)
                new_val = val;
            else
                key_set = keys(F.count_dict);
                cols = setdiff(F.current_set, e);
                ind1 = F.current_set ~= e;
                ind2 = logical([ind1, 1]);
                for i = 1:length(key_set)
                    key_string = char(key_set(i));
                    current_key = str2num(key_string);
                    if length(current_key) == n
                        ind = ind1;
                    else
                        ind = ind2;
                    end
                    rows = F.index_dict(key_string);
                    count = F.count_dict(key_string);
                    new_key = num2str(current_key(ind));
                    try
                        new_count_dict(new_key) = new_count_dict(new_key) + count;
                        new_index_dict(new_key) = [new_index_dict(new_key), rows];
                    catch
                        new_count_dict(new_key) = count;
                        new_index_dict(new_key) = rows;
                    end
                end
                
                F.count_dict = new_count_dict;
                F.index_dict = new_index_dict;
               
                ent = 0;
                joint = 0;
                key_set = keys(F.count_dict);
                m = n-1;

                for j = 1:length(key_set)
                    key = char(key_set(j));
                    key_array = str2num(key);
                    if length(key_array) == m
                        prob = F.count_dict(key)/F.n_rows;
                        ent = ent - (prob) * log2(prob);
                    else
                        prob = F.count_dict(key)/F.n_rows;
                        joint = joint - (prob) * log2(prob);
                    end
                end
                new_val = -1*(ent - joint) - F.classEntropy;
                
                F.current_set = cols;
                F.current_val = new_val;
                
            end
        end
   end
end