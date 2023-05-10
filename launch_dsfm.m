function [discrete_sol, discrete_obj, local_minimality, H, V, info] = launch_dsfm(dataset,lambda,method,param,seed) 
    % set seed
    rng(seed)
    % load data and corresponding difference of submodular function
    [H, V] = load_data_diffsub(dataset, lambda, seed);
    
    function [solution, bundle] = submin_algo(H, bundle)
        if param.warm_start
            if nargin < 3 
                [solution,~,~,~,~,bundle] = minimize_submodular_FW_minnormpoint_restart(H, V, param.inner_maxiter, 1, param.inner_gap);
            else
                [solution,~,~,~,~,bundle] = minimize_submodular_FW_minnormpoint_restart(H, V, param.inner_maxiter, 1, param.inner_gap, bundle);
            end
        else 
            solution = minimize_submodular_FW_minnormpoint(H, V, param.inner_maxiter, 1, param.inner_tol);
            bundle = [];
        end
    end

    switch method
        case "pgm" % ignore that objective is not submodular
            n = length(V);
            [discrete_sol, discrete_obj, info.gaps, info.itertime] = prox_lovasz_hypercube(H.H, V, zeros(n,1), 0, param.maxiter, param.gap, 1, zeros(n,1), Inf, true);
        case "mnp" % ignore that objective is not submodular
             [discrete_sol,~,~,discrete_obj,~,~,~,info.itertime] = minimize_submodular_FW_minnormpoint(H.H, V, param.maxiter, 1, param.gap);
        case "greedy" % max -H(S), ignore that objective is not supermodular
            invH = SetFctLinComb({H.H}, -1);
            [discrete_sol, discrete_obj, info.itertime] = Greedy_USM(invH,V, 1);
            discrete_obj = -discrete_obj;
        case "regdc" 
            prox = @(F, param_F, z, lambda, x_init) prox_lovasz_hypercube(F, V, z, lambda, param.inner_maxiter, param.inner_tol, 1, x_init, F(V), false);
            [discrete_sol, info.continuous_sol, discrete_obj, info.continuous_obj, info.gaps, info.perm_choices, info.itertime] = regDC(H, V, param.rho, param.perms, prox, param.maxiter, param.tol, param.warm_start, false, false, 0, param.localmin);
        case "regadc" 
            prox = @(F, param_F, z, lambda, x_init) prox_lovasz_hypercube(F, V, z, lambda, param.inner_maxiter, param.inner_tol, 1, x_init, F(V), false);
            [discrete_sol, info.continuous_sol, discrete_obj, info.continuous_obj, info.gaps, info.perm_choices, info.itertime] = regDC(H, V, param.rho, param.perms, prox, param.maxiter, param.tol, param.warm_start, false, true, param.q, param.localmin);
        case "regdcRound"
            prox = @(F, param_F, z, lambda, x_init) prox_lovasz_hypercube(F, V, z, lambda, param.inner_maxiter, param.inner_tol, 1, x_init, F(V), false);
            [discrete_sol, info.continuous_sol, discrete_obj, info.continuous_obj, info.gaps, info.perm_choices, info.itertime] = regDC(H, V, param.rho, param.perms, prox, param.maxiter, param.tol, param.warm_start, true, false, 0, param.localmin);
        case "regadcRound"
            prox = @(F, param_F, z, lambda, x_init) prox_lovasz_hypercube(F, V, z, lambda, param.inner_maxiter, param.inner_tol, 1, x_init, F(V), false);
            [discrete_sol, info.continuous_sol, discrete_obj, info.continuous_obj, info.gaps, info.perm_choices, info.itertime] = regDC(H, V, param.rho, param.perms, prox, param.maxiter, param.tol, param.warm_start, true, true, param.q, param.localmin);
        case "regcdc" 
            prox = @(F, param_F, z, lambda, x_init) prox_lovasz_hypercube(F, V, z, lambda, param.prox_maxiter, param.prox_tol, 1, x_init, F(V), false);
            [discrete_sol, info.continuous_sol, discrete_obj, info.continuous_obj, info.gaps, info.prox_max_gaps, info.perm_choices, info.itertime, info.FW_niters] = regCompleteDC(H, V, param.rho, param.perms, prox, param.maxiter,  param.tol, param.fw_maxiter, param.fw_gap, param.warm_start, param.ls, false, param.localmin);
        case "regcdcRound" 
            prox = @(F, param_F, z, lambda, x_init) prox_lovasz_hypercube(F, V, z, lambda, param.prox_maxiter, param.prox_tol, 1, x_init, F(V), false);
            [discrete_sol, info.continuous_sol, discrete_obj, info.continuous_obj, info.gaps, info.prox_max_gaps, info.perm_choices, info.itertime, info.FW_niters] = regCompleteDC(H, V, param.rho, param.perms, prox, param.maxiter,  param.tol, param.fw_maxiter, param.fw_gap, param.warm_start, param.ls, true, param.localmin);
        case "subsup" % equivalent to DC
             [discrete_sol, discrete_obj, info.perm_choices, info.itertime] = SubSup(H, V, param.perms, @submin_algo, param.maxiter, param.localmin);
        case "supsub" 
             % solve unconstrained sub max problem using algo of Buchbinder et al, 2012
             submax_algo = @(H) Greedy_USM(H, V, 1);
             [discrete_sol, discrete_obj, info.itertime] = SupSub(H, V, param.modFcts, submax_algo, param.maxiter, param.localmin);
        case "modmod" 
             [discrete_sol, discrete_obj, info.itertime] = ModMod(H, V, param.modFcts, param.perms, param.maxiter, param.localmin); 
        case "bruteforce"
            [discrete_sol, discrete_obj] = BruteForceMin(H.H, V);
            info = [];
    end
    
    % Check if solution is a discrete local min
    % set current sol of H to solution, for faster computation of marginals
    local_minimality = check_localmin(H, discrete_sol, V);
end