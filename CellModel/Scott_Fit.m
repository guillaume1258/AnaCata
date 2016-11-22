%%
%% Compute score
% Initialize matrices for lambda and the ratio RNA/total_prot (mean &
% std)
Scott.lambda_mean   = zeros(5 , 6);
Scott.RNA_Prot_mean = zeros(5 , 6);
Scott.lambda_std    = zeros(5 , 6);
Scott.RNA_Prot_std  = zeros(5 , 6);

% Mean  lambda
Scott.lambda_mean(: , 1)   = [0.40 , 0.33 , 0.24 , 0.19 , 0.12];
Scott.lambda_mean(: , 2)   = [0.57 , 0.50 , 0.39 , 0.30 , 0.23];
Scott.lambda_mean(: , 3)   = [0.71 , 0.57 , 0.38 , 0.23 , 0.14];
Scott.lambda_mean(: , 4)   = [1.00 , 0.87 , 0.67 , 0.43 , 0.28];
Scott.lambda_mean(: , 5)   = [1.31 , 0.90 , 0.46 , 0.20 , 0.11];
Scott.lambda_mean(: , 6)   = [1.58 , 1.18 , 0.89 , 0.31 , 0.13];

% lambda std
Scott.lambda_std(: , 1)    = [0.03 , 0.01 , 0.01 , 0.03 , 0.01];
Scott.lambda_std(: , 2)    = [0.02 , 0.02 , 0.04 , 0.02 , 0.03];
Scott.lambda_std(: , 3)    = [0.03 , 0.03 , 0.01 , 0.04 , 0.01];
Scott.lambda_std(: , 4)    = [0.05 , 0.05 , 0.04 , 0.06 , 0.05];
Scott.lambda_std(: , 5)    = [0.07 , 0.13 , 0.01 , 0.04 , 0.03];
Scott.lambda_std(: , 6)    = [0.15 , 0.10 , 0.08 , 0.12 , 0.02];

% RNA/prot mean
Scott.RNA_Prot_mean(: , 1) = [0.177 , 0.291 , 0.375 , 0.414 , 0.631];
Scott.RNA_Prot_mean(: , 2) = [0.230 , 0.303 , 0.371 , 0.400 , 0.496];
Scott.RNA_Prot_mean(: , 3) = [0.224 , 0.313 , 0.435 , 0.473 , 0.524];
Scott.RNA_Prot_mean(: , 4) = [0.287 , 0.340 , 0.374 , 0.471 , 0.577];
Scott.RNA_Prot_mean(: , 5) = [0.414 , 0.476 , 0.618 , 0.715 , 0.785];
Scott.RNA_Prot_mean(: , 6) = [0.466 , 0.500 , 0.584 , 0.691 , 0.769];

% RNA/prot std
Scott.RNA_Prot_std(: , 1)  = [0.006 , 0.014 , 0.015 , 0.028 , 0.092];
Scott.RNA_Prot_std(: , 2)  = [0.014 , 0.009 , 0.009 , 0.072 , 0.051];
Scott.RNA_Prot_std(: , 3)  = [0.029 , 0.037 , 0.045 , 0.042 , 0.079];
Scott.RNA_Prot_std(: , 4)  = [0.009 , 0.012 , 0.015 , 0.028 , 0.013];
Scott.RNA_Prot_std(: , 5)  = [0.058 , 0.063 , 0.081 , 0.065 , 0.100];
Scott.RNA_Prot_std(: , 6)  = [0.033 , 0.037 , 0.050 , 0.114 , 0.029];

%% Convert the data in the unit of measures in WS model

p = 0.76; % Conversion parameter to get phi_r from RNA/total_prot

%
Scott.alpha_mean = Scott.lambda_mean ./ 60 * log(2);
Scott.alpha_std  = Scott.lambda_std  ./ 60 * log(2);

%
Scott.f_r_mean = Scott.RNA_Prot_mean .* p;
Scott.f_r_std  = Scott.RNA_Prot_std  .* p;

Scott.alpha = mean(Scott.f_r_mean(: , :)) / mean(Scott.alpha_mean(: , :)) .* Scott.alpha_mean;
v.lambda    = v.lambda * 60 / log(2);