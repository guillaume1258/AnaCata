function [dX , flag , new_data] = ode( t , X , data)
%function [dX] = ode( t , X , k , dxdt)
% Integration occurs here

k = data.k;

%% Map vector X indices into specie names

Initial_derivatives;

S = Ind2Species(X);

%% Auxiliary functions of the state

A = Auxiliary_functions(k , S);

%% Write the CORE ODEs (corresponding to WS gene expression reactions)

dxdt = Core(S , k , A , dxdt);

%% Metabolism

        % sugar import
        dxdt.s_i = + dxdt.s_i + A.mu_import - A.mu_catalysis;
        % sugar metabolism
        dxdt.a   = + dxdt.a +  A.mu_catalysis * k.n_s;

%% Map the object dxdt.(names) to the dxdt vector indices

dX = Ind2Derivatives(dxdt) - X * A.lambda;

flag = 0;

new_data = [];

end

