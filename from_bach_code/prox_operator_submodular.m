function [w,bundle] = prox_operator_submodular(z,lambda,F,param_F,gap,bundle)
% minimize lambda * f(w) + .5 * || w - z ||^2
% using the min norm point algorithm
% allows random restarts (important in proximal methods)

if nargin<5, gap = 1e-8; end
param_Fprox.F = F;
param_Fprox.param_F = param_F;
param_Fprox.n = param_F.n;
param_Fprox.lambda = lambda;
param_Fprox.z = z;

if nargin == 6
    % modify bundle to take into account the shift in z and lambda
    bundle.x = bundle.x * lambda - z;
    bundle.X = bundle.X * lambda - repmat(z,1,size(bundle.X,2));
    [w,t1,t2,t3,t4,bundle] = minimize_submodular_FW_minnormpoint_restart(@submodular_fct_prox,param_Fprox,param_F.n*2,0,gap,bundle);
else
    [w,t1,t2,t3,t4,bundle] = minimize_submodular_FW_minnormpoint_restart(@submodular_fct_prox,param_Fprox,param_F.n*2,0,gap);
end
    figure, plot(t3)
% modify bundle to take into account the shift in z and lambda
    bundle.x =  ( z+ bundle.x ) / lambda ;
    bundle.X = ( bundle.X + repmat(z,1,size(bundle.X,2)) ) / lambda;

  