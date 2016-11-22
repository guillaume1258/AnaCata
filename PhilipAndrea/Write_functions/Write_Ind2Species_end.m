function [  ] = Write_Ind2Species_end( species )
% Write the function that gives the value to the species names from the
% corresponding index in the vector X

%% Number of species

N_s = size(species , 2);

%% Write the function Ind2Species

New = fopen('Ind2Species_end.m' , 'w');

fprintf(New , 'function [S] = Ind2Species_end(X) \n');
fprintf(New , '%% gives the value to the species names from the corresponding index in the vector X \n');

for i = 1 : N_s
    
    s = species{i};
    
    fprintf(New , 'S.%s(: , 1) = X(: , %d); \n' , s , i);
    
end

fprintf(New , 'end');


end
