%% Load and plot results
close all
clear all
addpath(genpath(pwd));
expdir = "speech-lbd1.0e+00-d268453-53048204/"; %"speech-lbd1.0e+00-d268453-53048204/"; %"mushroom-lbd1.0e-04-2dc9f47-56258106/"; %"speech-lbd1.0e+00-173be10-14674247/"; "speech-lbd1.0e+00-7fc4444-48639887/"; %"mushroom-lbd1.0e-04-2dc9f47-56258106/"; %"speech-lbd1.0e+00-d268453-53048204/";  %; %mushroom-lbd1.0e-03-2dc9f47-56258081 
loaddir ="results/" + expdir;  
figsdir = "figs/" + expdir;
save_fig = true;
if save_fig
    mkdir(figsdir);
end
files = dir(fullfile(loaddir,'*.mat'));
%% load results

methods_run = [];

for i=1:length(files)
    load(loaddir + files(i).name);
    ind_run = seed - 42 + 1;
    total_results.(method_name)(ind_run) = results.(method_name);
    if ~ismember(method_name, ["mnp", "greedy", "pgm"]) %
        % remove last iteration 
        total_results.(method_name)(ind_run).discrete_obj = total_results.(method_name)(ind_run).discrete_obj(1:end-1);
        total_results.(method_name)(ind_run).itertime = total_results.(method_name)(ind_run).itertime(1:end-1);
        if isfield(total_results.(method_name)(ind_run), 'continuous_obj')
            total_results.(method_name)(ind_run).continuous_obj = total_results.(method_name)(ind_run).continuous_obj(1:end-1);
            total_results.(method_name)(ind_run).gaps = total_results.(method_name)(ind_run).gaps(1:end-1);
        end
    elseif method_name=="mnp"
         total_results.(method_name)(ind_run).discrete_obj = total_results.(method_name)(ind_run).discrete_obj';
         total_results.(method_name)(ind_run).itertime = total_results.(method_name)(ind_run).itertime';
    end
    if i==1
        methods_run = [method];
    else
        methods_run = union(methods_run, method);
    end

    if isfield(param.(method_name), 'fw_maxiter')
        fw_maxiter = param.(method_name).fw_maxiter;
        total_results.(method_name)(ind_run).FW_niters = total_results.(method_name)(ind_run).FW_niters(1:end-1);
        total_results.(method_name)(ind_run).prox_max_gaps = total_results.(method_name)(ind_run).prox_max_gaps(1:end-1);
    end
end
%methods_run = setdiff(methods_run, ["mnp", "modmod", "supsub", "greedy"]); %temp ignore problematic cases
%% get min objectives
best_rho = struct();
total_results.min_continuous  = zeros(1, nruns);
total_results.min_discrete = zeros(1, nruns);
for method=methods_run
    method_name = strrep(method, '.', '');
    rho_str = regexp(method,'\d+\.?\d*','Match');
    method_class = erase(method, rho_str);
    for ind_run = 1:length(total_results.(method_name))
        total_results.(method_name)(ind_run).min_discrete = min(total_results.(method_name)(ind_run).discrete_obj);
        total_results.min_discrete(ind_run) = min([total_results.min_discrete(ind_run), total_results.(method_name)(ind_run).min_discrete]);
        if ismember(method_class, ["regdc", "regdcRound", "regadc", "regadcRound", "regcdc", "regcdcRound", "pgm"]) %
                total_results.(method_name)(ind_run).min_continuous = min(total_results.(method_name)(ind_run).continuous_obj);
                total_results.min_continuous (ind_run) = min([total_results.min_continuous (ind_run), total_results.(method_name)(ind_run).min_continuous]);
        end
    end
    if ismember(method_class, ["regdc", "regdcRound", "regadc", "regadcRound", "regcdc", "regcdcRound"])
        [val, ind] =  min(mean(padcat(total_results.(method_name).discrete_obj),2, 'omitnan'));
        if ~isfield(best_rho, method_class) ||  val < best_rho.(method_class).val || (val == best_rho.(method_class).val && ind < best_rho.(method_class).ind)
            best_rho.(method_class).rho = rho_str;
            best_rho.(method_class).val = val;
            best_rho.(method_class).ind = ind;
        end
    end
  
end
total_results.min_continuous
total_results.min_discrete
best_rho
%% Check running times
% fprintf("Running times: \n")
% for method=methods_run
%     method_name = strrep(method, '.', '');
%     for ind_run = 1:length(total_results.(method_name))
%         if isfield(total_results.(method_name), 'time')
%             fprintf("%s: %s \n", method, total_results.(method_name).time)
%         elseif ~isempty(total_results.(method_name)(ind_run).itertime) %isfield(total_results.(method_name)(ind_run), 'itertime')
%             fprintf("%s: %s \n", method, total_results.(method_name)(ind_run).itertime(end))
%         end
%     end 
% end


%% Check perm choices
fprintf("Perm choices: \n")
for method=methods_run
    method_name = strrep(method, '.', '');
    for ind_run = 1:length(total_results.(method_name))
        if isfield(total_results.(method_name)(ind_run), 'perm_choices')
            fprintf("%s : \n", method_name)
            total_results.(method_name)(ind_run).perm_choices'
        end
    end
end
%% Check local minimality

fprintf("Local minimality: \n") % should be negative for local min sol

if ~exist('H') || ~ exist('V')
    rng(seed)
    [H, V] = load_data_diffsub(dataset, lambda);
end

for method=methods_run
    method_name = strrep(method, '.', '');
    for ind_run = 1:length(total_results.(method_name))
        if isfield(total_results.(method_name)(ind_run), 'local_minimality')
            fprintf("%s: %e \n", method, total_results.(method_name)(ind_run).local_minimality)
        else % get local minimality for results where we did not save it
            solution = total_results.(method_name).solution;
            [H_sol, H.H] = H.H(solution);
            % double check that H_sol == last discrete obj
            total_results.(method_name).discrete_obj(end)
            H_sol
            min_H_local = 0;
            for i= setdiff(V, solution)
                min_H_local = min(min_H_local, add(H.H, solution, i));
            end
            for i= solution(:)' 
                min_H_local = min(min_H_local, rmv(H.H, solution, i));
            end
            local_minimality = H_sol - min_H_local;
            fprintf("%s: %e \n", method, local_minimality)
        end
    end
end

%% plot results

labels_map = containers.Map(["regdc", "regadc", "regdcRound", "regadcRound", "regcdc", "regcdcRound", "subsup","supsub", "modmod", "mnp", "greedy", "pgm"], ["DCA", "ADCA", "DCAR", "ADCAR", "CDCA", "CDCAR", "SubSup", "SupSub", "ModMod", "MNP", "Greedy", "PGM"]);
include_zero = true;
eps = 1e-6; % to avoid log(0)
delta = 3;
if dataset == "speech"
    delta_time = 2.5e3;
    xlimit = 100;
else
    delta_time = 1e4; %2.5e3; %1e4
    xlimit = 150;
end
%time_intervals = 0:2.5e3:2e4;
%% plot discrete or continous objective vs iteration or time (only include best rho for DCA and CDCA)
discrete = true;
only_baselines = false;
only_zerorho = true;
if discrete
    classes = ["regdc","regadc", "regdcRound", "regadcRound", "regcdc", "regcdcRound", "pgm", "subsup","supsub", "modmod", "mnp", "greedy"];
else
    classes = ["regdc","regadc", "regdcRound", "regadcRound", "regcdc", "regcdcRound", "pgm"];
end
colors = distinguishable_colors(length(classes)-2); 
colors = [colors(1, :); colors(1, :);colors(2, :); colors(2, :); colors(3:end, :)];
lines = ["-x", "--x", "-*", "--*", "-o", "-s", "-p", "-h", "-d", "-^", "-+", "-v"];
figure%('Position', [0, 0, 600, 400])
hold all
labels = [];
xaxis = "time (sec)"; %  "iterations", "time (sec)"

if only_baselines
    classes_ind = [7, 9:length(classes)];
    delta_time = 10;
elseif only_zerorho
    classes_ind = [1,3,8];
else
    classes_ind = 1:length(classes);
end
for i = classes_ind
    method_class = classes(i);
    if ismember(method_class, ["regdc", "regdcRound", "regadc", "regadcRound", "regcdc", "regcdcRound"]) && isfield(best_rho, method_class)
        if only_zerorho
            rho_str = "0";
        else
            rho_str = best_rho.(method_class).rho;
        end
        method = method_class+ rho_str;
        method_name = strrep(method, '.', '');
        method_label = labels_map(method_class)+ " $\rho =$ " + rho_str;
    else
        rho_str = "";
        method = method_class;
        method_name = method_class;
        method_label = labels_map(method_class);
    end

    if ~ismember(method, methods_run)
        continue
    end
    [x, y] = get_xy(method_name, method_class, discrete, total_results, include_zero, xaxis, outer_maxiter, eps);
    if xaxis ==  "time (sec)"
        indices = sample_intervals(x, 0:delta_time:x(end)+delta_time);   
    else
        indices = sample_intervals(x, 0:delta:x(end)+delta); % to handle regcdc methods, equivalent to 1:delta:length(x) for other methods
    end
    labels = [labels, method_label];
    plot(x, mean(y, 2, 'omitnan'), lines(i), 'color', colors(i, :),'linewidth',2, 'markerindices', indices,'markersize',10); 
end

for i = classes_ind
    method_class = classes(i);
    if ismember(method_class, ["regdc", "regdcRound","regadc", "regadcRound", "regcdc", "regcdcRound"]) &&  isfield(best_rho, method_class)
        if only_zerorho
            rho_str = "0";
        else
            rho_str = best_rho.(method_class).rho;
        end
        method = method_class+ rho_str;
        method_name = strrep(method, '.', '');
    else
        method = method_class;
        method_name = method_class;
    end
    if ~ismember(method, methods_run)
        continue
    end
    
    [x, y] = get_xy(method_name, method_class, discrete, total_results, include_zero, xaxis, outer_maxiter, eps);

    if xaxis ==  "time (sec)"
        indices = sample_intervals(x, 0:delta_time:x(end)+delta_time);   
    else
        indices = sample_intervals(x, 0:delta:x(end)+delta);
    end
    x = x(indices);
    y = y(indices, :);
    errorbar(x, mean(y, 2, 'omitnan'), min(std(y, 0, 2, 'omitnan'), mean(y, 2, 'omitnan')-eps), '.', 'color', colors(i, :), 'linewidth', 1, 'Capsize', 0);

end

set(gca,'fontsize',20,'YScale', 'log', 'XScale', 'linear')
xlabel(xaxis)
if discrete
    ylabel('$F(X^k) - \min(F)$', 'Interpreter','latex') % we use F instead of H to match notation in paper.
else
    ylabel('$f_L(x^k) - \min(f_L)$', 'Interpreter','latex')
end
if only_baselines 
    xlim([0, xlimit])
else
    l = legend(labels, 'Location','northeast');
    set(l,'Interpreter','latex', 'fontsize',15)
    axis tight
end

%set(gcf, 'PaperUnits', 'centimeters','PaperPosition', [0 0 22 15]); %
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

if save_fig
    if discrete
        if xaxis == "iterations"
            fig_name = figsdir + sprintf("/discrete-obj-iter-bestrho");
        else
            fig_name = figsdir + sprintf("/discrete-obj-time-bestrho");
        end
        if only_baselines 
                fig_name = fig_name + "-baselines";
        elseif only_zerorho
            fig_name = fig_name + "-zerorho";
        end
    else
        if xaxis == "iterations"
            fig_name = figsdir + sprintf("/cont-obj-iter-bestrho");
        else
            fig_name = figsdir + sprintf("/cont-obj-time-bestrho"); 
        end
        if only_baselines 
            fig_name = fig_name + "-baselines";
        elseif only_zerorho
            fig_name = fig_name + "-zerorho";
        end
    end
    print(gcf,'-dpdf','-r150',fig_name); %'-bestfit'
else
    plotbrowser
end

%% plot discrete or continuous objective of our methods with all rhos
classes = ["regdc", "regdcRound", "regadc", "regadcRound", "regcdc", "regcdcRound"];
rhos = [0, 0.001, 0.01, 0.1, 1, 10];
nrhos = length(rhos);
colors = distinguishable_colors(nrhos); 
lines = ["-x","-*", "-o", "-s", "-p", "-h"];
xaxis =  "iterations"; %  "iterations", "time (sec)"
discrete = true;
figure('Position', 0.8*get(0, 'Screensize'))
index = reshape(1:length(classes),3,2).';
if ~islogical(discrete) && discrete == "gap"
    eps = 1e-6;
    include_zero = false;
end
for i = 1:length(classes)
    labels = [];
    subplot(2,3,index(i));
    hold all
    method_class = classes(i);
    for ind_rho = 1:nrhos
        rho_str = num2str(rhos(ind_rho));
        method = method_class+ rho_str;
        method_name = strrep(method, '.', '');
        if ~ismember(method, methods_run)
            continue
        end
        
        [x, y] = get_xy(method_name, method_class, discrete, total_results, include_zero, xaxis, outer_maxiter, eps);
  
        labels = [labels,  labels_map(method_class)+ " $\rho =$ " + rho_str];
        if xaxis == "time (sec)"
            indices = sample_intervals(x, 0:delta_time:x(end)+delta_time);   
        else
            indices = sample_intervals(x, 0:delta:x(end)+delta);
        end
        plot(x, mean(y, 2, 'omitnan'), lines(ind_rho), 'color', colors(ind_rho, :),'linewidth',2, 'markerindices', indices,'markersize',10); 
    end
    
    for ind_rho = 1:nrhos
        rho_str = num2str(rhos(ind_rho));
        method = method_class+ rho_str;
        method_name = strrep(method, '.', '');
        if ~ismember(method, methods_run)
            continue
        end
        
        [x, y] = get_xy(method_name, method_class, discrete, total_results, include_zero, xaxis, outer_maxiter, eps);
        
        if xaxis == "time (sec)"
            indices = sample_intervals(x, 0:delta_time:x(end)+delta_time);   
        else
            indices = sample_intervals(x, 0:delta:x(end)+delta);
        end
        x = x(indices);
        y = y(indices, :);
        errorbar(x, mean(y, 2, 'omitnan'), min(std(y, 0, 2, 'omitnan'), mean(y, 2, 'omitnan')-eps), '.', 'color', colors(ind_rho, :), 'linewidth', 1, 'Capsize', 0);
    end
    l = legend(labels, 'Location','northeast');
    set(l,'Interpreter','latex', 'fontsize',15)
    set(gca,'fontsize',20,'YScale', 'log', 'XScale', 'linear')
    xlabel(xaxis)
    if ~islogical(discrete) && discrete=="gap"
        ylabel('gap', 'Interpreter','latex')
    elseif discrete
        ylabel('$F(X^k) - \min(F)$', 'Interpreter','latex') % we use F instead of H to match notation in paper.
    else
        ylabel('$f_L(x^k) - \min(f_L)$', 'Interpreter','latex')
    end
    axis tight
end

%set(gcf, 'PaperUnits', 'centimeters','PaperPosition', [0 0 35 30]); %
set(gcf,'Units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

if save_fig
     if ~islogical(discrete) && discrete=="gap"
         fig_name = figsdir + sprintf("/gap-obj-iter-allrhos");
     elseif discrete 
        if xaxis == "iterations"
            fig_name = figsdir + sprintf("/discrete-obj-iter-allrhos");
        else
            fig_name = figsdir + sprintf("/discrete-obj-time-allrhos");
        end
    else
        if xaxis == "iterations"
            fig_name = figsdir + sprintf("/cont-obj-iter-allrhos");
        else
            fig_name = figsdir + sprintf("/cont-obj-time-allrhos");
        end
    end
    print(gcf,'-dpdf','-r100',fig_name); %'-bestfit'
    %saveas(gcf,fig_name, '-dpdf','-r150');
else
    plotbrowser
end


%% functions for plotting

function [x, y] = get_xy(method_name, method_class, discrete, total_results, include_zero, xaxis, outer_maxiter, eps)

    if ~islogical(discrete) && discrete=="gap"
        min_obj = zeros(1,3); 
        if ismember(method_class, ["regcdc", "regcdcRound"])
            y = padcat(total_results.(method_name).prox_max_gaps);
        else
            y = padcat(total_results.(method_name).gaps);
        end
    elseif discrete
        min_obj = total_results.min_discrete;
        y = padcat(total_results.(method_name).discrete_obj);
    else
        min_obj = total_results.min_continuous;
        y = padcat(total_results.(method_name).continuous_obj);
    end
    nseeds_run = size(y,2); % temporary fix for some methods not running for all seeds
    y = y - min_obj(1:nseeds_run) + eps;

    if xaxis ==  "time (sec)"
            x =  mean(padcat(total_results.(method_name).itertime), 2, 'omitnan');
    else
        x = (1:size(y,1))';
    end
    
    if include_zero
            y = [zeros(1, nseeds_run) - min_obj(1:nseeds_run) + eps; y];
            x = [0; x];
    end
    
    if ismember(method_class, ["regcdc", "regcdcRound"]) && (xaxis == "iterations") && islogical(discrete)
            x = cumsum(mean(padcat(total_results.(method_name).FW_niters), 2, 'omitnan'));
            if include_zero
                x = [0; x];
            end 

    elseif ismember(method_class,["mnp", "modmod", "supsub", "greedy", "pgm", "bruteforce"]) && (xaxis == "iterations")
            x = (1:outer_maxiter)';
            if discrete
                y = [total_results.(method_name).min_discrete] .* ones(outer_maxiter, 1) -  min_obj(1:nseeds_run) + eps;
            else
                y = [total_results.(method_name).min_continuous] .* ones(outer_maxiter, 1) -  min_obj(1:nseeds_run) + eps;
            end
            if include_zero
                 y = [y(1,:); y];
                 x = [0; x];
            end
    end
end

function indices = sample_intervals(data, intervals)
    indices = [];
    for i=1:length(intervals)-1
        indices = [indices, find(data >= intervals(i) & data <= intervals(i+1), 1)];
    end
end


%% Older plots 
% 
% %% plot discrete objective
% 
% colors = distinguishable_colors(length(methods_run)); 
% figure
% hold all
% eps = 1e-6; % to avoid log(0)
% disp("Methods achieving best discrete obj: ")
% % for i = 1:length(methods_run) % ugly fix to remove error bars from legend
% %     plot(0,0,'color', colors(i, :), 'linewidth',2);
% % end
% for i = 1:length(methods_run)
%     method = methods_run(i);
%     method_name = strrep(method, '.', '');
%     rho_str = regexp(method,'\d+\.?\d*','Match');
%     method_class = erase(method, rho_str);
%     if ismember(method_class, ["subsup","regdc", "regdcRound", "regadc", "regadcRound"]) % we don't plot value at beginning (zero for all) for clarity 
%             y = [zeros(1, nruns); total_results.(method_name).discrete_obj] - total_results.min_discrete + eps;
%             errorbar(0:size(y,1)-1, mean(y, 2), std(y, 0, 2),'color', colors(i, :), 'linewidth',2, 'Capsize', 0);
%             % plot(1e-6+total_results.(method_name).discrete_obj - total_results.min_discrete,'color', colors(i, :), 'linewidth',2);
%     elseif ismember(method_class, ["regcdc", "regcdcRound"])
%             y = [total_results.(method_name).discrete_obj];
%             temp = size(y,2); % temporary fix to some methods not running for all seeds
%             y = [zeros(1, temp); y] - total_results.min_discrete(1:temp) + eps; 
%             if isfield(total_results.(method_name), 'FW_niters')
%                 errorbar(cumsum(mean([zeros(1, temp); total_results.(method_name).FW_niters], 2)), mean(y, 2), std(y, 0, 2),'color', colors(i, :),'linewidth',2, 'Capsize', 0);   
%             else
%                 niter = length(total_results.(method_name)(1).discrete_obj);
%                 errorbar(fw_maxiter:fw_maxiter:min(maxiter,niter*fw_maxiter), mean(y, 2), std(y, 0, 2),'color', colors(i, :),'linewidth',2, 'Capsize', 0);   
%             end
%     elseif ismember(method_class,["mnp", "modmod", "supsub", "greedy", "bruteforce"])
%             y = [zeros(1, nruns); [total_results.(method_name).min_discrete] .* ones(outer_maxiter/5, nruns)] -  total_results.min_discrete + eps;
%             errorbar(0:5:outer_maxiter, mean(y, 2), std(y, 0, 2), 'color', colors(i, :),'linewidth',2, 'Capsize', 0); 
%             % plot(1e-6+ min(total_results.(method_name).discrete_obj) * ones(outer_maxiter, 1) - total_results.min_discrete, 'color', colors(i, :),'linewidth',2); 
%     end
%     %fprintf("%s: %f\n", method, min(total_results.(method_name).discrete_obj))
% %     if total_results.(method_name).min_discrete <= total_results.min_discrete + 1e-5
% %         disp(method)
% %     end
% end
% %[lgd, icons, plots, txt] = legend(methods_run);
% l = legend(methods_run, 'Location','NorthEastOutside');
% set(gca,'fontsize',25,'YScale', 'log', 'XScale', 'linear')
% xlabel('iterations')
% %ylabel('F(S)') 
% ylabel('$F(S) - F^\star$', 'Interpreter','latex') 
% %axis tight
% xlim([1,inf])
% set(l,'Interpreter','latex')
%  
% if save_fig
%     fig_name = loaddir + sprintf("/discrete-obj-outiter-%s", dataset);
%     print(gcf,'-dpdf','-r150',fig_name);
% else
%     plotbrowser
% end
% 
% %% plot continuous objective
% figure
% hold all
% %fprintf("Best continuous objective value achieved by: \n")  
% disp("Methods achieving best continuous obj: ")
% for i = 1:length(methods_run)
%     method = methods_run(i);
%     method_name = strrep(method, '.', '');
%     rho_str = regexp(method,'\d+\.?\d*','Match');
%     method_class = erase(method, rho_str);
%     if ismember(method_class, ["regdc", "regdcRound", "regadc", "regadcRound"])
%             plot(total_results.(method_name).continuous_obj ,'color', colors(i, :), 'linewidth',2); 
%             if total_results.(method_name).min_continuous <= min_continuous + 1e-5
%                 disp(method)
%             end
%     elseif ismember(method_class, ["regcdc", "regcdcRound"])
%             niter = length(total_results.(method_name).continuous_obj);
%             if isfield(total_results.(method_name), 'FW_niters')
%                 plot(cumsum(total_results.(method_name).FW_niters), total_results.(method_name).continuous_obj, 'color', colors(i, :),'linewidth',2); 
%             else
%                 plot(fw_maxiter:fw_maxiter:min(maxiter,niter*fw_maxiter), total_results.(method_name).continuous_obj, 'color', colors(i, :),'linewidth',2);   
%             end
%     end
%     %fprintf("%s: %f\n", method,  min(total_results.(method_name).continuous_obj))
% 
% end
% legend(setdiff(methods_run, ["mnp", "modmod", "subsup", "supsub", "greedy", "bruteforce"], 'stable')) 
% set(gca,'fontsize',20,'YScale', 'log', 'XScale', 'linear')
% xlabel('iterations')
% ylabel('h_L(x)')
% %ylabel('h_L(x) - min h_L')
% axis tight
% 
% % save fig
% if save_fig
%     fig_name = loaddir + sprintf("/continuous-obj-outiter-%s", dataset);
%     print(gcf,'-dpdf','-r150',fig_name);
% end
% plotbrowser
% 
% %%
% % plot gaps
% figure
% hold all
% for i = 1:length(methods_run)
%     method = methods_run(i);
%     method_name = strrep(method, '.', '');
%     rho_str = regexp(method,'\d+\.?\d*','Match');
%     method_class = erase(method, rho_str);
%     if ismember(method_class, ["regdc", "regdcRound", "regadc", "regadcRound"])
%             plot(total_results.(method_name).gaps,'color', colors(i, :), 'linewidth',2); 
%     elseif ismember(method_class, ["regcdc", "regcdcRound"])
%             niter = length(total_results.(method_name).gaps);
%             if isfield(total_results.(method_name), 'FW_niters')
%                 plot(cumsum(total_results.(method_name).FW_niters), total_results.(method_name).gaps, 'color', colors(i, :),'linewidth',2);
%             else
%                 plot(fw_maxiter:fw_maxiter:min(maxiter,niter*fw_maxiter), total_results.(method_name).gaps, 'color', colors(i, :),'linewidth',2);   
%             end
%     end
%     %fprintf("%s: %f\n", method,  min(total_results.(method_name).continuous_obj))
% 
% end
% legend(setdiff(methods_run, ["mnp", "modmod", "subsup", "supsub", "greedy", "bruteforce"], 'stable')) 
% set(gca,'fontsize',20,'YScale', 'linear', 'XScale', 'linear')
% xlabel('iterations')
% ylabel('gap(x)')
% axis tight
% 
% % save fig
% if save_fig
%     fig_name = loaddir + sprintf("/gaps-outiter-%s", dataset);
%     print(gcf,'-dpdf','-r150',fig_name);
% end
% plotbrowser
% %%
% fprintf("CDCA max prox gaps: \n")
% for method=intersect(methods_run, [regcdc_variants, regcdcRound_variants])
%     method_name = strrep(method, '.', '');
%     if isfield(total_results.(method_name), 'prox_max_gaps')
%         fprintf("%s: ", method)
%         %total_results.(method_name).prox_max_gaps'
%         total_results.(method_name).FW_niters'
%     end
% end
% %% plot discrete objective vs time
% figure
% hold all
% for i = 1:length(methods_run)
%     method = methods_run(i);
%     method_name = strrep(method, '.', '');
%     plot(total_results.(method_name).itertime, total_results.(method_name).discrete_obj, 'color', colors(i, :),'linewidth',2); 
% end
% legend(methods_run)
% set(gca,'fontsize',20,'YScale', 'log', 'XScale', 'linear')
% xlabel('time (sec)')
% ylabel('F(S)')
% %xlim([1,inf])
% if save_fig
%     fig_name = filename + sprintf("/discrete-obj-outiter-%s", dataset);
%     print(gcf,'-dpdf','-r150',fig_name);
% end
% plotbrowser
% 
