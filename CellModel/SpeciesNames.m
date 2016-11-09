function [ species ] = SpeciesNames( Gene_names , Type )
%SPECIESNAMES Summary of this function goes here
%   Detailed explanation goes here

%% Number of genes

N_g = size(Gene_names , 1);

%% Number of types

N_t = size(Type , 1);

%% Index gene related species

for i = 1 : N_t
    
    TypeIs = Type{i};
    
    for j = 1 : N_g
        
        gene = Gene_names{j};
        
        species{(i - 1) * N_g + j} = sprintf('%s_%s' , TypeIs , gene);
        
    end
    
end

%% Index non gene related species

species{N_t * N_g + 1} = 's_i';
species{N_t * N_g + 2} = 'a';
%species{N_t * N_g + 3} = 'e_DnaA_a';

end

