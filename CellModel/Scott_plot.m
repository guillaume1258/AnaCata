function [ ] = Scott_plot( v )

% WS model data
WS_Model = importdata('WS_Scott_Data.mat');

% Make Scott figure

% Mean  lambda
Exp_Data.lambda_mean(: , 1)   = [0.40 , 0.33 , 0.24 , 0.19 , 0.12];
Exp_Data.lambda_mean(: , 2)   = [0.57 , 0.50 , 0.39 , 0.30 , 0.23];
Exp_Data.lambda_mean(: , 3)   = [0.71 , 0.57 , 0.38 , 0.23 , 0.14];
Exp_Data.lambda_mean(: , 4)   = [1.00 , 0.87 , 0.67 , 0.43 , 0.28];
Exp_Data.lambda_mean(: , 5)   = [1.31 , 0.90 , 0.46 , 0.20 , 0.11];
Exp_Data.lambda_mean(: , 6)   = [1.58 , 1.18 , 0.89 , 0.31 , 0.13];

% lambda std
Exp_Data.lambda_std(: , 1)    = [0.03 , 0.01 , 0.01 , 0.03 , 0.01];
Exp_Data.lambda_std(: , 2)    = [0.02 , 0.02 , 0.04 , 0.02 , 0.03];
Exp_Data.lambda_std(: , 3)    = [0.03 , 0.03 , 0.01 , 0.04 , 0.01];
Exp_Data.lambda_std(: , 4)    = [0.05 , 0.05 , 0.04 , 0.06 , 0.05];
Exp_Data.lambda_std(: , 5)    = [0.07 , 0.13 , 0.01 , 0.04 , 0.03];
Exp_Data.lambda_std(: , 6)    = [0.15 , 0.10 , 0.08 , 0.12 , 0.02];

% RNA/prot mean
Exp_Data.RNA_Prot_mean(: , 1) = [0.177 , 0.291 , 0.375 , 0.414 , 0.631];
Exp_Data.RNA_Prot_mean(: , 2) = [0.230 , 0.303 , 0.371 , 0.400 , 0.496];
Exp_Data.RNA_Prot_mean(: , 3) = [0.224 , 0.313 , 0.435 , 0.473 , 0.524];
Exp_Data.RNA_Prot_mean(: , 4) = [0.287 , 0.340 , 0.374 , 0.471 , 0.577];
Exp_Data.RNA_Prot_mean(: , 5) = [0.414 , 0.476 , 0.618 , 0.715 , 0.785];
Exp_Data.RNA_Prot_mean(: , 6) = [0.466 , 0.500 , 0.584 , 0.691 , 0.769];

% RNA/prot std
Exp_Data.RNA_Prot_std(: , 1)  = [0.006 , 0.014 , 0.015 , 0.028 , 0.092];
Exp_Data.RNA_Prot_std(: , 2)  = [0.014 , 0.009 , 0.009 , 0.072 , 0.051];
Exp_Data.RNA_Prot_std(: , 3)  = [0.029 , 0.037 , 0.045 , 0.042 , 0.079];
Exp_Data.RNA_Prot_std(: , 4)  = [0.009 , 0.012 , 0.015 , 0.028 , 0.013];
Exp_Data.RNA_Prot_std(: , 5)  = [0.058 , 0.063 , 0.081 , 0.065 , 0.100];
Exp_Data.RNA_Prot_std(: , 6)  = [0.033 , 0.037 , 0.050 , 0.114 , 0.029];

% Convert the data in the unit of measures in WS model
p = 0.76; % Conversion parameter to get phi_r from RNA/total_prot

%
Exp_Data.phi_r_mean = Exp_Data.RNA_Prot_mean .* p;
Exp_Data.phi_r_std  = Exp_Data.RNA_Prot_std  .* p;

figure(1)
hold all

cmap = hsv(6);
for i = 1 : 6
    %plot(v.lambda(: , i)  , v.phi_r(: , i) , 'LineWidth' , 3);
    WS_model_plot = plot(WS_Model.lambda(: , i) * 60 / log(2), WS_Model.phi_r(: , i)  ,  'x-.' , 'Color' , cmap(i , :) , 'LineWidth' , 2);  
    model_plot = plot(v.lambda(: , i) , v.phi_r(: , i)  ,  'o-' , 'Color' , cmap(i , :) , 'LineWidth' , 2);    
    Data_plot = plot(Exp_Data.lambda_mean(: , i) ./ log(2) , Exp_Data.phi_r_mean(: , i) , 'o' , 'Color' , cmap(i , :), 'LineWidth' , 3);
end

ylim([0 0.6])

x_label = xlabel('Growth-rate [1/hour]');
y_label = ylabel('Ribosomal mass fraction \Phi_r');
title_s = title('w_r = 2000');
%title_s = title('w_r = 930');

set(y_label ,'FontSize',16);
set(x_label ,'FontSize',16);
set(title_s ,'FontSize',16);

set(gca , 'FontSize' , 14);

legend([WS_model_plot , model_plot , Data_plot] , 'WS Model' , 'Same n_x' , 'Data')

Save_Name = 'Figures/Scott_laws_Reduced_wr_2000.eps';

print(gcf , Save_Name , '-dpsc2')
end

