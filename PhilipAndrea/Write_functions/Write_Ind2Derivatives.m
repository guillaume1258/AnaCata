function [  ] = Write_Ind2Derivatives( species )
% Write the function that gives a value to the Derivatives 
% from the corresponding index in the vector dxdt

%% Number of species

N_s = size(species , 2);

%% Write the function Ind2Derivatives

New = fopen('Ind2Derivatives.m' , 'w');

fprintf(New , 'function [dX] = Ind2Derivatives(dxdt) \n');
fprintf(New , '%% gives the value to the species names from the corresponding index in the vector X \n');


for i = 1 : N_s
    
    s = species{i};
    
    fprintf(New , 'dX(%d , 1) = dxdt.%s; \n' , i , s);
    
end




end


