function f = submodular_fct_prox(A,param_F)
% affine transform of a submodular function
% used within proximal operators
f = param_F.F(A,param_F.param_F) * param_F.lambda - sum( param_F.z(A) );