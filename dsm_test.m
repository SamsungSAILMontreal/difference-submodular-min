function dsm_test(nruns, dataset, lambda, task_id, job_id, resultsdir)
% runs dsm for one method and one seed
    disp([class(task_id), class(job_id)]);
    seeds = 42:42+nruns;
    addpath(genpath(pwd));
    [~,git_hash_string] = system('git rev-parse --short HEAD');
    dirname = sprintf(resultsdir + "/%s-lbd%1.1e-%s-%s",dataset, lambda,erase(git_hash_string,newline), num2str(job_id));
    mkdir(dirname);
    %%
    regdc_variants = [];
    regadc_variants = [];
    regdcRound_variants = [];
    regadcRound_variants = [];
    regcdc_variants = [];
    regcdcRound_variants = [];
    for rho = [0, 0.001, 0.01, 0.1, 1, 10]
        regdc_variants = [regdc_variants, "regdc" + num2str(rho)];
        regadc_variants = [regadc_variants, "regadc" + num2str(rho)];
        regdcRound_variants = [regdcRound_variants, "regdcRound" + num2str(rho)];
        regadcRound_variants = [regadcRound_variants, "regadcRound" + num2str(rho)];
        regcdc_variants = [regcdc_variants, "regcdc" + num2str(rho)];
        regcdcRound_variants = [regcdcRound_variants, "regcdcRound" + num2str(rho)];
    end
    methods = [regdc_variants, regadc_variants, regdcRound_variants, regadcRound_variants, regcdc_variants, regcdcRound_variants, "modmod", "subsup", "supsub", "mnp", "greedy", "pgm"];
    maxiter = 3e4; 
    inner_tol = 1e-6;
    outer_tol = 1e-6;
    inner_maxiter = 1e3;
    outer_maxiter = ceil(maxiter/inner_maxiter);
    % time of 1 mnp iter = time of greedy algo (sorting + n EO) + proj on
    % affine hull (n^3? to be checked)
    % time of 1 subsup iter = time of greedy algo + submin prob (computation for 
    % perms done in parallel)
    % time of 1 inner iter of subsup (submin prob) = time of 1 mnp iteration
    % time of 1 regdc iter = time of greedy algo + prox prob (computation for 
    % perms done in parallel)
    % time of 1 inner iter of regdc (prox prob)= time of greedy algo 
    % time of 1 regcdc iter = FW + prox prob
    % time of 1 inner iter of regcdc (FW)= time of greedy algo + prox prob
%%
    method = methods(mod(task_id , length(methods))+1);
    seed = seeds(fix(task_id/length(methods))+1);
    fprintf('Running %s on %s dataset with lambda=%f and seed=%d: \n',method, dataset,lambda,seed)
    method_name = strrep(method, '.', '');
    rho_str = regexp(method,'\d+\.?\d*','Match');
    method_class = erase(method, rho_str);
    param.(method_name) = struct;
    if ismember(method, ["subsup", regdc_variants, regadc_variants, regdcRound_variants, regadcRound_variants])
        param.(method_name).inner_maxiter = inner_maxiter;
        param.(method_name).inner_tol = inner_tol;
        param.(method_name).maxiter = outer_maxiter;
        param.(method_name).perms = [1,2,3];
        param.(method_name).tol = outer_tol;
        param.(method_name).localmin = true;
        if ismember(method, [regadc_variants, regadcRound_variants])
            param.(method_name).q = 5;
        end
        if method=="subsup" %TODO: check why we get negative duality gap in mnp with warm start m
            param.(method_name).warm_start = false;
        else
            param.(method_name).warm_start = true;
            param.(method_name).rho = str2num(rho_str);
        end
        
     elseif ismember(method, [regcdc_variants, regcdcRound_variants])
        param.(method_name).fw_maxiter = 10;
        param.(method_name).inner_maxiter = inner_maxiter; 
        param.(method_name).maxiter = outer_maxiter; %ceil(maxiter/(param.(method_name).fw_maxiter*inner_maxiter));
        param.(method_name).perms = [1,2,3];
        param.(method_name).rho = str2num(rho_str);
        param.(method_name).fw_gap = inner_tol;
        param.(method_name).inner_tol = inner_tol;
        param.(method_name).tol = outer_tol;
        param.(method_name).ls = false;
        param.(method_name).warm_start = true; 
        param.(method_name).localmin = true;
        
    elseif ismember(method, ["mnp", "modmod", "supsub", "pgm"])
        param.(method_name).maxiter = maxiter;
        param.(method_name).gap = -inf; %no gap stopping criteria %outer_tol;
        param.(method_name).perms = [1,2,3];
        param.(method_name).modFcts = [1,2];
        param.(method_name).localmin = true;
    end
    
    %tic;
    [discrete_sol, discrete_obj, local_minimality, H, V, info] = launch_dsfm(dataset,lambda,method_class,param.(method_name),seed);
    %time = toc;
    
    results.(method_name).discrete_sol = discrete_sol;
    results.(method_name).discrete_obj = discrete_obj;
    results.(method_name).local_minimality = local_minimality;
    for fn = fieldnames(info)'
        results.(method_name).(fn{1}) = info.(fn{1});
    end
    %results.(method_name).time = time;
    %% Save results
    filename = dirname + sprintf("/%d-%s-%s", seed, method_name, datetime('now'));
    fprintf("saving results to %s", filename)
    save(filename);
end