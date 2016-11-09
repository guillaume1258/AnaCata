function [  ] = Write_Initial_Conditions( species )
% Write the function that define initial conditions

%% Number of species

N_s = size(species , 2);

%% Write the function that define initial conditions

New = fopen('Initial_Conditions.m' , 'w');

%% Define the script

fprintf(New , '%% Set the initial conditions of the model \n \n');

for i = 1 : N_s
    
    s = species{i};
    
    fprintf(New , '%s = 0; \n' , s);
    
end


fprintf(New , 'e_r = 1000; \n');
fprintf(New , 'e_t = 1000; \n');
fprintf(New , 's_i = 0; \n');
fprintf(New , 'a   = 10000; \n \n');


for i = 1 : N_s
    
    s = species{i};
    fprintf(New , 'X0(%d) = %s; \n' , i , s);    
    
end


end

