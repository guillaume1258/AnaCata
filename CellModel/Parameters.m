k = struct();                                                                    

alpha = 1; %3 * 1e4;

k.beta  = 1.63;

%% Masses

k.n = struct(); 

        % cell proteic mass
        k.M     = 10^8;
    
        % individual protein masses
        k.n.r     = 7459 * (1 + k.beta);         % length of ribosome [aa / molecs]
        k.n.t     = 300;          % length of non ribosomal T proteins [aa/ molecs]
        k.n.m     = 300;          % length of non ribosomal M proteins [aa/ molecs]
        k.n.p     = 300;          % length of non ribosomal P proteins [aa/ molecs]
        k.n.q     = 300;          % length of non ribosomal Q proteins [aa/ molecs]        
        
%% Expression

k.w     = struct();
k.theta = struct();
k.K     = struct();

        % transcriptional parameters
        k.w.r     = 31;          % max ribosome transcription rate [molecs / min / cell]
        k.w.t     = 3.4;         % max enzyme T transcription rate [molecs / min /cell]
        k.w.m     = 3.4;         % max enzyme M transcription rate [molecs / min /cell]
        k.w.p     = 0;            % max enzyme P transcription rate [molecs / min /cell]
        k.w.q     = 840.93;       % max q transcription rate [molecs / min /cell]
        
        k.K.q     = 100000;                     % q auto-inhibition threshold [molecs / cell]
        k.alpha_q = 4;                          % q auto-inhibition hill coefficient [none] 
        
        k.theta.r   = 4570 * alpha;       % ribosome transcription threshold [molecs / cell]
        k.theta.t   = 240   * alpha;       % t protein transcription threshold [molecs / cell]
        k.theta.m   = 240   * alpha;       % m protein transcription threshold [molecs / cell]
        k.theta.p   = 240   * alpha;       % p protein transcription threshold [molecs / cell]
        k.theta.q   = 240   * alpha;       % q protein transcription threshold [molecs / cell]
        
        k.rif       = 0;                            % rifampicin concentration [microM]
        k.K.rif     = 3;                         % Rifampicin action threshold [microM]
                
        k.d         = 0.1;          % [r,t,m,p,q,...] mRNA-degradation rate [min^-1]
        
    % translational parameters
        k.gamma_max = 1260;         % max translational elongation rate [aa / min / molecs]
        k.K.p       = 7.8;          % scaling translational elongation threshold [aa.cell/min]
        k.K.gamma = k.gamma_max / k.K.p * alpha;
        
%% competition

        % mRNA-ribosome binding parameters
        k.b       = 1;            % r mRNA-ribosome binding rate [cell / min / molecs ]
        
        % mRNA-ribosome unbinding parameters
        k.u       = 1;            % r mRNA-ribosome unbinding rate [min^-1]
        
        % mRNA-ribosome binding to chloramphenicol parameters
        k.b_cm    = 0.00599;      % chloramphenicol-binding rate [(min.microM)^-1] 
        k.cm      = 0;            % quantity of chloramphenicol
        
%% Metabolism

k.v = struct();

        k.v.t = 726;                   % max. nutrient import rate [min^-1]
        k.K.t = 1000;                  % nutrient import threshold [molecs]
        k.v.m = 5800;                  % max enzymatic rate [min^-1]
        k.K.m = 1000;                  % enzymatic threshold [molecs / cell]
        
        k.n_s = 0.5;                   % nutrient efficiency [none]

% %% Toggle-switch
% 
%         % Transcription
%         k.w.a     = 80; 
%         k.w.b     = 10;
%         k.theta.a = 10;
%         k.theta.b = 10;
%         k.K.a     = 100;
%         k.K.b     = 10;
%         k.alpha_a = - 1;
%         k.alpha_b = - 1;
%         
%         % Translation
%         k.n.a     = 300;
%         k.n.b     = 300;
%         
%         % Competition
%         k.k.b.a   = 1;
%         k.k.b.b   = 1;
%         k.k.u.a   = 1;
%         k.k.u.b   = 1;
%         
%         % Toggle inhibition
%         k.I_a     = 1;
%         k.I_b     = 1;
%         k.k.b.I_a = 1;
%         k.k.b.I_b = 1;
%         k.k.u.I_a = 1;
%         k.k.u.I_b = 1;
        
%% environment

        k.s_0 = 10^4;
        
%% Replication

        k.b_DnaA = 1;
        k.u_DnaA = 1;