
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

N_Lim  = 4;
N_Cond = 6;

k.sigma_i = 1;

for i = 1 : length(v.n_s)
    
    for j = 1 : 1
        
        Parameters;
        k.n_s = v.n_s(i);
        
        %X0(1:26) = X_Init(i , j).X';
%         
%         TxScaling = v.TxScaling(j);
%         
%         k.w.r = k.w.r * TxScaling;
%         k.w.t = k.w.t * TxScaling;
%         k.w.m = k.w.m * TxScaling;
%         k.w.p = k.w.p * TxScaling;
%         k.w.q = k.w.q * TxScaling;
%         
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

save('SimulationsData/CurrentRawData.mat' , 'SteadyX')