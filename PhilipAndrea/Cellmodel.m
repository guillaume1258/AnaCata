
%% Clear work memory

close all
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

Gene_names = {'r' ; 't' ; 'm' ; 'q' ; 'p'};


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
X_Init = importdata('X0_vectors/Initial_Condition_Hwa.mat');
%% Loop over range of n_s and [cm]

v.n_s    = [0.106 , 0.116 , 0.218 , 0.366 , 0.564 , 1.123]; 
v.cm     = [0 , 0 , 2 , 4 , 6 , 8];
v.wp     = [0 , 10.4 , 35.6 , 72.1 , 213 , 395];
v.wp_cAA = [0 , 17.6 , 40 , 60.7 , 130 , 210 , 299];
v.TxScaling = logspace(- 3 , 0 , 20); 

v.n_s = linspace(0.1 , 1.5 , 10);

Nut_Lim.sigma_i     = [1.3054 , 1.2774 , 1.3380 , 1.4795 , 1.5997 , 2.1006];
CHL_Lim.sigma_i     = [1.4795 , 1.4795 , 1.4757 , 1.4891 , 1.5012 , 1.5266];
OE_Lim.sigma_i      = [1.5011 , 1.5011 , 1.5583 , 1.6736 , 1.9213 , 2.0918];
OE_CasA_Lim.sigma_i = [1.7120 , 1.7885 , 1.7892 , 2.0461 , 2.1072 , 2.1054];

N_Lim  = 4;
N_Cond = 6;

k.sigma_i = 1;

for i = 1 : length(v.n_s)
    
    for j = 1 : 1%length(v.cm)
        
        Parameters;
        k.n_s = v.n_s(i);
        
        %X0(1:26) = X_Init(i , j).X';
        %% Define various event counters
%         
%         if i == 1
%             
%             k.n_s     = v.n_s(j);
%             k.cm      = 0;
%             k.w.p     = 0;
%             k.sigma_i = Nut_Lim.sigma_i(j);
%             
%         elseif i == 2
%             
%             k.cm      = v.cm(j);
%             k.n_s     = v.n_s(4);
%             k.w.p     = 0;
%             k.sigma_i = CHL_Lim.sigma_i(j);
%             
%         elseif i == 3
%             
%             k.cm      = 0;
%             k.n_s     = v.n_s(4);
%             k.w.p     = v.wp(j);
%             k.sigma_i = OE_Lim.sigma_i(j);
%             
%         elseif i == 4
%             k.cm      = 0;
%             k.n_s     = v.n_s(5);
%             k.w.p     = v.wp_cAA(j + 1);
%             k.sigma_i = OE_CasA_Lim.sigma_i(j);
%             
%             
%         end
        
        TxScaling = v.TxScaling(j);
        TxScaling = 7;
        
        k.w.r = k.w.r * TxScaling;
        k.w.t = k.w.t * TxScaling;
        k.w.m = k.w.m * TxScaling;
        k.w.p = k.w.p * TxScaling;
        k.w.q = k.w.q * TxScaling;
        
        tic;
        k.StartTime = clock;
        %% Integration time
        
        t0    = 0;
        tf    = 1e6;
        tspan = [t0 tf];
        
        %% CVode
        data = struct();
        data.k = k;
        
        options = CVodeSetOptions('UserData', data ,...
            'RelTol',1.e-6,...
            'AbsTol',1e-9*ones(size(X0)),...
            'LinearSolver','Dense',...
            'InitialStep',1e-3,...
            'MaxStep',1000,...
            'MaxNumSteps',1e5);
        CVodeInit(@ode, 'BDF', 'Newton', t0, X0' , options);
        [status,t,X] = CVode([tf],'Normal');
        
        % Record the steady state
        SteadyX(i , j).X       = X;
        
        % Get the number of i protein for each condition
        a_number(i , j)        = X(end);
        
        toc;
        
    end
    
end

