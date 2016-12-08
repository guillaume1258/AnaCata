
k = struct();

% rRNA multiplicator for ribosome's mass and cost.
k.beta  = 1.63;

%% Masses

k.n = struct();

% cell proteic mass
k.M     = 10^8;
k.M_0   = 10^8;
k.M_ref = 10^8;

% individual protein masses
k.n.r     = 7459;         % length of ribosome [aa / molecs]
k.n.t     = 300;          % length of non ribosomal T proteins [aa/ molecs]
k.n.m     = 300;          % length of non ribosomal M proteins [aa/ molecs]
k.n.p     = 300;          % length of non ribosomal P proteins [aa/ molecs]
k.n.q     = 300;          % length of non ribosomal Q proteins [aa/ molecs]

%% Expression

k.w     = struct();
k.theta = struct();
k.K     = struct();

z = 1;
% transcriptional parameters
k.w.r     = 1   * z;         % max ribosome transcription rate [molecs / min / cell]
k.w.t     = 0.4 * z;         % max enzyme T transcription rate [molecs / min /cell]
k.w.m     = 0.4 * z;         % max enzyme M transcription rate [molecs / min /cell]
k.w.p     = 0;                 % max enzyme P transcription rate [molecs / min /cell]
k.w.q     = 18.98 * z;         % max q transcription rate [molecs / min /cell]

k.K.q     = 167030;            % q auto-inhibition threshold [molecs / cell]
k.alpha_q = 4;                 % q auto-inhibition hill coefficient [none]

k.theta.r   = 5180;       % ribosome transcription threshold [molecs / cell]
k.theta.t   = 0.41;       % t protein transcription threshold [molecs / cell]
k.theta.m   = 0.41;       % m protein transcription threshold [molecs / cell]
k.theta.p   = 0.41;       % p protein transcription threshold [molecs / cell]
k.theta.q   = 0.41;       % q protein transcription threshold [molecs / cell]

k.d         = 0.1;          % [r,t,m,p,q,...] mRNA-degradation rate [min^-1]

% translational parameters
k.gamma_max = 12600;         % max translational elongation rate [aa / min / molecs]
k.K.p       = 40.8;          % scaling translational elongation threshold [aa.cell/min]
k.K.gamma = k.gamma_max / k.K.p;

%% competition

% mRNA-ribosome binding parameters
k.b       = 1;            % r mRNA-ribosome binding rate [cell / min / molecs ]

% mRNA-ribosome unbinding parameters
k.u       = 1;            % r mRNA-ribosome unbinding rate [min^-1]

% mRNA-ribosome binding to chloramphenicol parameters
k.b_cm    = 0.0033;      % chloramphenicol-binding rate [(min.microM)^-1]
k.cm      = 0;            % quantity of chloramphenicol

%% Metabolism

k.v = struct();

k.v.t = 726;                   % max. nutrient import rate [min^-1]
k.K.t = 1000;                  % nutrient import threshold [molecs]
k.v.m = 5800;                  % max enzymatic rate [min^-1]
k.K.m = 1000;                  % enzymatic threshold [molecs / cell]

k.n_s = 0.5;                   % nutrient efficiency [none]

%% PMF

k.PMF = -170; % Membrane potential, in mV

%% environment

k.s_0 = 10^4;
