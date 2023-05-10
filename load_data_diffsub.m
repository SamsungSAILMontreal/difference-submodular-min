function [H, V] = load_data_diffsub(dataset, lambda, seed)
% load data and define H difference of submodular functions F - G s.t
% H.H(S, param_H) = H.F(S, param_H.F)- H.G(S, param_H.G)
% param_H should have field param_H.n = size of ground set

switch dataset
    case "speech"
        W = read_bilmes('uniform-sent.800');
        
        % choose a small subset of data for faster testing
        % W = W(1:50, :);
        % randomize order
        rp = randperm(size(W,1));
        W = W(rp,:);
        
        V = 1:size(W,1);
        
        alpha = 0.5;
        F = SetFctLinComb({CoverFct(W, alpha)}, lambda);
        
        n = size(W,1);
        m = 10;
        % partition ground set into m subsets
        grp_edges = round(linspace(1,n+1,m+1));
        % sample random non-negative weights
        c = abs(randn(n, 1));
        G = SubPartitionFct(c, grp_edges);
        
        H = diff_sub_fct(F, G);
%         % compute beta of G and its lower bound
%         beta_lb = 1;
%         for i=1:m
%             c_i = c(grp_edges(i):grp_edges(i+1)-1);
%             beta_lb = min(beta_lb, 0.5*(min(c_i)/sum(c_i))^0.5);
%         end
%         fprintf("beta lower bound= %f \n", beta_lb);
%         
%         beta = 1;
%         [GV, G] = G(V); 
%         for i = V
%             beta = min(beta, ( GV - rmv(G, V, i))/G([i]));
%         end
%         fprintf("beta= %f \n", beta);
        
    case "mushroom"
        W = readmatrix(sprintf('datasets/mushroom_train_%d.csv', seed-42+1));
        %test with small subset of data
        %W = W(:,1:15);
        %separate class variable from data
        C = W(:,1);
        %remove class from data
        W(:,1) = [];
        V = 1:size(W,2);
        F = ConditionalEntropyReg(W,C, lambda);
        G = Entropy(W); %child of SetFct
        H = MutualInformationReg(W,C,lambda); %child of SetFct
        H = diff_sub_fct(F, G, H);
        
    case "adult"
        W = readmatrix('datasets/adult_train.csv');
        %test with small subset of data
        %W = W(1:100,[1,4:10]);
        %separate class variable from data
        C = W(:,1);
        %remove class from data
        W(:,1) = [];
        V = 1:size(W,2);
        F = ConditionalEntropyReg(W,C, lambda);
        G = Entropy(W); %child of SetFct
        H = MutualInformationReg(W,C,lambda); %child of SetFct
        H = diff_sub_fct(F, G, H);
        
end
