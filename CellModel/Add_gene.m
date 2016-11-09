function [ Gene_names , k ] = Add_gene( Gene_names , k , new_gene , theta , w )
% Add the transcription parameters for a gene and add its component to the
% Genes_names vector

%% Size of Gene_names

N_g = size(Gene_names , 1);

Gene_names{end + 1} = new_gene;
k.theta.(new_gene)  = theta;
k.w.(new_gene)      = w;
k.n.(new_gene)      = 300;

end

