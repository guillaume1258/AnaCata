function [  ] = Write_Initial_derivatives( species )
% Write the Initial_derivatives function

%% Number of species

N_s = size(species , 2);

%% Write the script Initial_derivatives

New = fopen('Initial_derivatives.m' , 'w');

fprintf(New , '%% Set default value of the derivatives to 0 \n');

for i = 1 : N_s
    
    s = species{i};
    
    fprintf(New , 'dxdt.%s = 0; \n' , s);
    
end

end

