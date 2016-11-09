function [ ] = Write_CORE( Gene_names )
% Automatically write the reaction of the core model for a set of genes


%% Number of genes

N_s = size(Gene_names , 1);


%% Create the ODEs corresponding to the core for those genes

New = fopen('Core.m' , 'w');


%% Define the function

fprintf(New , 'function [dxdt] = Core(S , k , A , dxdt) \n');
fprintf(New , '%% ODEs corresponding to the CORE gene equations \n \n');

%% Compute differential equations

for i = 1 : N_s
    
    gene = Gene_names{i};
    
    fprintf(New , '\n \n');
    
    % -> m_x @omega_x
    if gene == 'q'
    fprintf(New , 'dxdt.m_q = + dxdt.m_q + Hill(1 , k.theta.q , k.w.q , S.a) * Hill(k.alpha_q , S.e_q , 1 , k.K.q); \n');
    else
    fprintf(New , 'dxdt.m_%s = + dxdt.m_%s + Hill(1 , k.theta.%s , k.w.%s , S.a); \n \n' , gene , gene , gene , gene);
    end
    
    % m_x -> @d
    fprintf(New , 'dxdt.m_%s = + dxdt.m_%s - k.d * S.m_%s; \n \n' , gene ,gene , gene);
    
    
    % n_x * a + rm_x -> e_r + m_x + e_x @Nu_x
    fprintf(New , 'v1_%s      = + Nu(S.rm_%s , 300 , A.gamma); \n' , gene , gene);
    fprintf(New , 'dxdt.m_%s  = + dxdt.m_%s + v1_%s; \n' , gene , gene , gene);
    fprintf(New , 'dxdt.e_r   = + dxdt.e_r  + v1_%s; \n' , gene); 
    fprintf(New , 'dxdt.e_%s  = + dxdt.e_%s + 300 / k.n.%s * v1_%s; \n' , gene , gene , gene , gene);
    fprintf(New , 'dxdt.rm_%s = + dxdt.rm_%s - v1_%s; \n' , gene , gene , gene);
    fprintf(New , 'dxdt.a     = + dxdt.a - v1_%s * 300; \n \n' , gene);
       
    
    % e_r + m_x <-> rm_x @k.b , k.u
    fprintf(New , 'v2_%s = - S.m_%s * S.e_r * k.b + k.u * S.rm_%s; \n' , gene , gene , gene);
    fprintf(New , 'dxdt.m_%s = + dxdt.m_%s + v2_%s; \n' , gene , gene , gene);
    fprintf(New , 'dxdt.e_r = + dxdt.e_r + v2_%s; \n' , gene);
    fprintf(New , 'dxdt.rm_%s = + dxdt.rm_%s - v2_%s; \n \n' , gene , gene , gene);
    
    % rm_x -> zm_x @[cm]*k_cm
    fprintf(New , 'v3_%s = + S.rm_%s * k.cm * k.b_cm; \n' , gene , gene);
    fprintf(New , 'dxdt.rm_%s = + dxdt.rm_%s - v3_%s; \n' , gene , gene , gene);
    fprintf(New , 'dxdt.zm_%s = + dxdt.zm_%s + v3_%s; \n \n' , gene , gene , gene);
         
end

fprintf(New , '\n');

fprintf(New , 'end');

end

