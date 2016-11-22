function [  ] = Write_Auxiliary_functions( Gene_names )
% Write the auxiliary functions of the model

%% Number of genes

N_g = size(Gene_names , 1);

%% Define the function

New = fopen('Auxiliary_functions.m' , 'w');

fprintf(New , 'function [A] = Auxiliary_functions(k , S) \n');
fprintf(New , '%% Compute the various functions of the state \n \n');

fprintf(New , 'A.mu_import    = mu(S.e_t , 1 , k.K.t , k.v.t , k.s_0) ; \n');   % Total transport rate
fprintf(New , 'A.mu_catalysis = mu(S.e_m , 1 , k.K.m , k.v.m , S.s_i); \n');    % Total metabolic rate

fprintf(New , 'A.gamma        = Hill(1 , k.K.gamma , k.gamma_max , S.a); \n');  % Translation elongation rate

fprintf(New , 'A.rate_tot = (');

for i = 1 : N_g
    
    fprintf(New , '+ S.rm_%s ' , Gene_names{i});
    
end

fprintf(New , '+ k.beta * S.rm_r) * A.gamma; \n');

fprintf(New , 'A.lambda       = A.rate_tot / k.M; \n');

fprintf(New , 'end');

end

