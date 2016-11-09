%% Clear work memory

clear all


startup_STB;

%% Add paths

addpath('Generic_functions/');
addpath('Write_functions/');

%% Load parameters

Parameters;

%% Define type of species

Type = {'m' ; 'rm' ; 'zm' ; 'e'};

%% CORE genes

Gene_names = {'r' ; 't' ; 'm' ; 'q'};

%% Add new genes

[Gene_names , k] = Add_gene(Gene_names , k , 'p' , 4.38 , 0); % Add gene p expression to the model


%% Names of the species

species = SpeciesNames(Gene_names , Type);

%% Write the functions that map the indices to the specie and derivative names

Write_Ind2Species(species);
Write_Ind2Species_end(species);
Write_Ind2Derivatives(species);

%% Write the auxiliary functions

Write_Auxiliary_functions(Gene_names);

%% Write the ODE file corresponding to CORE reactions

Write_CORE(Gene_names);

%% Write the function that set the default value of dxdt.(specie) to 0

Write_Initial_derivatives(species);

%% Write the function that set the initial conditions

Write_Initial_Conditions(species);

%% Load initial conditions

Initial_Conditions;
Initial_derivatives;

X0 = importdata('InitialCond.mat');

%% Loop over range of n_s and [cm]
v.cm     = [0 , 2 , 4 , 8 , 12];

v.n_s = logspace(log10(0.08) , log10(1) , 6);

X0 = X0';


for i = 1 : length(v.n_s)

    tic;
    k.n_s = v.n_s(i);
    
    for j = 1 : length(v.cm);
        k.cm  = v.cm(j);

        %% Integration time
        
        t0    = 0;
        tf    = 1e5;
        tspan = [t0 tf];
        
        %% Set options for the integration
  
%%
data = struct();
data.k = k;

        options = CVodeSetOptions('UserData', data,...
                         'RelTol',1.e-6,...
                         'AbsTol',1e-9*ones(size(X0)),...
                         'LinearSolver','Dense',...
                         'InitialStep',1e-3,...
                         'MaxStep',1000,...
                         'MaxNumSteps',1e5);
       CVodeInit(@ode, 'BDF', 'Newton', t0, X0 , options);
       [status,t,y] = CVode([tf],'Normal');

        
        %options  = odeset('NonNegative', [1:end] , 'RelTol' , 1e-6 , 'AbsTol' , 1e-9);              % Variables must be positive at all steps
        
        %% Run integrator
        
        %CVodeInit(@(t , X) ode(t , X , k , dxdt) , tspan , X0 , options);                   % Call ode solver
        [status,t,X]    = CVode(tf,'Normal');
        
        %% Load the X matrices into vectors with our human names
        
        S = Ind2Species_end(X');
        
        %% Compute observables
        
        obs = observables(S , k);
        
        v.phi_r(j , i)       = obs.r_fraction;
        v.lambda(j , i)      = obs.lambda;
        v.a(j , i)           = S.a(end);
        v.s_i(j , i)         = S.s_i(end);
        v.phi_t(j , i)       = obs.t_fraction(end);
        v.phi_q(j , i)       = obs.q_fraction(end);
        v.gamma(j , i)       = obs.gamma;
        v.mRNA_tot(j , i)    = S.m_r(end) + S.m_t(end) + S.m_m(end) + S.m_q(end) + S.m_p(end); 
        v.zm_rm_ratio(j , i) = (S.zm_r(end) + S.zm_t(end) + S.zm_m(end) + S.zm_q(end))...
                              /(S.rm_r(end) + S.rm_t(end) + S.rm_m(end) + S.rm_q(end) ...
                              + (S.zm_r(end) + S.zm_t(end) + S.zm_m(end) + S.zm_q(end)));
        
        v.C_D(j , i)         = 1 / (1.1880 / (2.0203 / v.a(j , i) + 1)) * 60;
        % on the side...
        %v.a_lam_ratio(j , i) = v.a ./ v.lambda;
        
        toc;
    end
    
    
end
%%
v.a_lam_ratio = v.a ./ v.lambda;

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

%%
Scott_plot(v)
