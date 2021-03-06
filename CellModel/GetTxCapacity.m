% Get the tx capacity in each condition

Parameters;

SteadyX = importdata('SimulationsData/CurrentRawData.mat');

v.n_s    = [0.106 , 0.116 , 0.218 , 0.366 , 0.564 , 1.123];
v.cm     = [0 , 0 , 2 , 4 , 6 , 8];
v.wp     = [0 , 10.4 , 35.6 , 72.1 , 213 , 395];
v.wp_cAA = [0 , 17.6 , 40 , 60.7 , 130 , 210 , 299];
v.TxScaling = logspace(- 3 , 0 , 20);


for i = 1 : 6
    
    k.n_s = v.n_s(i);
    
    % load Variable at stationary state and rename them
    X = SteadyX(i , j).X;
    S = Ind2Species(X);
    
    % Total tx activity for mRNA
    v.SumOmega(i , j) = k.n.r / k.n.q * k.w.r / (k.theta.r / S.a + 1) + ...
        k.w.t / (k.theta.t / S.a + 1) + ...
        k.w.m / (k.theta.m / S.a + 1) + ...
        k.w.p / (k.theta.p / S.a + 1) + ...
        k.w.q / (k.theta.q / S.a + 1);
    
    
    % sum tx with rRNA contribution: k.w.r / (k.theta.r / S.a + 1) * alpha
    alpha = 10;
    v.SumTx(i , j) = v.SumOmega(i , j) + k.w.r / (k.theta.r / S.a + 1) * alpha;
    
    % Tx fraction
    v.TxmRNAFraction(i , j) = v.SumOmega(i , j) / v.SumTx(i , j);
    
    % Max tx rate
    v.TxMax(i , j)    = k.w.r + k.w.t + k.w.m + k.w.p + k.w.q;
    
    % Total mRNA
    v.TotmRNA(i , j)  = S.m_r + S.m_t + S.m_m + S.m_p + S.m_q;
    
    % Total ribosome
    v.RTot(i , j)     = S.e_r + 10 * (S.rm_r + S.rm_t + S.rm_m + S.rm_p + S.rm_q + S.zm_r + S.zm_t + S.zm_m + S.zm_p + S.zm_q);
    
    % Non chl bound ribo
    v.c(i , j) = 10 * (S.rm_r + S.rm_m + S.rm_t + S.rm_p + S.rm_q) + S.e_r;
    
    % Tx capacity is the current tx activity divided by the max tx
    v.TxCapacity(i , j) = v.SumOmega(i , j) / v.TxMax(i , j);
    
    v.gamma(i , j)     = k.gamma_max / (k.K.gamma / S.a + 1);
    v.lambda(i , j)    = (S.rm_r * (1 + k.beta) + S.rm_t + S.rm_m + S.rm_q + S.rm_q) * v.gamma(i , j) / k.M_ref;
    
    % Ratio lambda tx capacity
    v.RTxLambda(i , j) = v.lambda(i , j) / v.SumOmega(i , j);
    
    % Free ribosomes
    v.FreeRibo(i , j)  = S.e_r;
    
    % Total proteic mass
    v.ProtTot(i , j)   = v.RTot(i , j) * k.n.r + S.e_t * k.n.t + S.e_m * k.n.m + S.e_p * k.n.p + S.e_q * k.n.q;
    
    % Total Mass
    v.Mass(i , j)      = v.ProtTot(i , j) + k.beta * k.n.r * v.RTot(i , j);
    
    % Mass fraction of proteins
    v.phi_r(i , j) = v.RTot(i , j) * k.n.r * (1 + k.beta) / (v.ProtTot(i , j) + v.RTot(i , j) * k.beta * k.n.r);
    v.phi_t(i , j) = S.e_t         * k.n.t / v.ProtTot(i , j);
    v.phi_m(i , j) = S.e_m         * k.n.m / v.ProtTot(i , j);
    v.phi_p(i , j) = S.e_p         * k.n.p / v.ProtTot(i , j);
    v.phi_q(i , j) = S.e_q         * k.n.q / v.ProtTot(i , j);
    
end

%% Plot fraction of free ribosomes as function of TxScaling factor (~Number of chromosomes)
figure(1)
hold all

for i = 1 : length(v.n_s)
    
    plot(log10(v.TxScaling) , v.lambda(i , :) / log(2) * 60 , 'LineWidth' , 3)
    
end

x_label = xlabel('log_{10}(Tx Max scaling)' , 'fontsize' , 18);
y_label = ylabel('\lambda [Doubling/hour]' , 'fontsize' , 18);

grid on

[hleg1, hobj1] = legend('n_s = 0.106' , 'n_s = 0.116' , 'n_s = 0.218' , 'n_s = 0.366' , 'n_s = 0.564' , 'n_s = 1.123');
set(hleg1, 'position' , [0.1 0.6 0.5 0.4])

set(gca , 'fontsize' , 18);

print(gcf , 'TxScaling_vs_lambda.eps' , '-dpsc2')

figure(2)
hold all

for i = 1 : length(v.n_s)
    
    plot(log10(v.TxScaling) , v.FreeRibo(i , :) ./ v.RTot(i , :) , 'LineWidth' , 3)
    
end


x_label = xlabel('log_{10}(Tx Max scaling)' , 'fontsize' , 18);
y_label = ylabel('Free Ribosome Fraction' , 'fontsize' , 18);

grid on

set(gca , 'fontsize' , 18);

print(gcf , 'TxScaling_vs_FreeRiboFraction.eps' , '-dpsc2')