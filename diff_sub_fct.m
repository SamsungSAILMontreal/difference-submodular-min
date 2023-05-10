function H = diff_sub_fct(F, G, FmG)
% Struct for difference of submodular functions FmG = F - G
% F, G, and FmG (optional) should be function objects of class sfo_fn  
H.F = F;
H.G = G;
if nargin<3
    H.H = SetFctLinComb({H.F,H.G},[1, -1]);
else
    H.H = FmG;
end
