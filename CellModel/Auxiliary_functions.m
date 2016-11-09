function [A] = Auxiliary_functions(k , S) 
% Compute the various functions of the state 
 
A.mu_import    = mu(S.e_t , 1 , k.K.t , k.v.t , k.s_0) ; 
A.mu_catalysis = mu(S.e_m , 1 , k.K.m , k.v.m , S.s_i); 
A.gamma        = Hill(1 , k.K.gamma , k.gamma_max , S.a); 
A.rate_tot = (+ S.rm_r + S.rm_t + S.rm_m + S.rm_q + S.rm_p ) * A.gamma; 
A.lambda       = A.rate_tot / k.M; 
end