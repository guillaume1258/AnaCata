function [ o ] = observables( S , k)
% Compute observables from the output variables of the model

%% observables

alpha = 1;

o = struct();

% total number of ribosomes
o.r_tot = 10 * (+ S.rm_r + S.rm_t + S.rm_m  + S.rm_q + S.rm_p ...
               + S.zm_r + S.zm_t + S.zm_m  + S.zm_q + S.zm_p) + S.e_r;

% total mass of each protein
o.ribosome_mass = ...
      k.n.r * (1 + k.beta) * o.r_tot;
  
o.t_mass = S.e_t * k.n.t;
o.m_mass = S.e_m * k.n.m;
o.q_mass = S.e_q * k.n.q;
o.p_mass = S.e_p * k.n.p;

% total mass
o.total_mass = ...
      o.ribosome_mass + ...
      o.t_mass + o.m_mass + o.q_mass + o.p_mass;

% rRNA mass
o.rRNA_mass = k.beta * o.ribosome_mass;
  
% mass fractions of r, r_free r, t, m, p and q  
o.r_fraction           = o.ribosome_mass(end) / o.total_mass(end);
o.r_free_fraction      = S.e_r * k.n.r ./ o.total_mass(end);
o.t_fraction           = S.e_t * k.n.t ./ o.total_mass(end);
o.m_fraction           = S.e_m * k.n.m ./ o.total_mass(end);
o.q_fraction           = S.e_q * k.n.q ./ o.total_mass(end);

% lambda
k.K.gamma    = k.gamma_max / k.K.p * alpha;
o.gamma      = Hill(1 , k.K.gamma , k.gamma_max , S.a(end));                                                         % Translation elongation rate
rate_tot     = (S.rm_r(end) * (1 + k.beta) + S.rm_t(end) + S.rm_m(end) + S.rm_q(end) + S.rm_p(end)) * o.gamma;  % Total polymerisation rate
o.lambda     = rate_tot / k.M;                                                                                       % Growth-rate

end

