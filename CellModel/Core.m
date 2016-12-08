function [dxdt] = Core(S , k , A , dxdt) 
% ODEs corresponding to the CORE gene equations 
 

 
dxdt.m_r = + dxdt.m_r + Hill(1 , k.theta.r , k.w.r , S.a); 
 
dxdt.m_r = + dxdt.m_r - k.d * S.m_r; 
 
v1_r      = + Nu(S.rm_r , k.n.r , A.gamma); 
dxdt.m_r  = + dxdt.m_r + v1_r; 
dxdt.e_r   = + dxdt.e_r  + 10 * v1_r; 
dxdt.e_r  = + dxdt.e_r + v1_r; 
dxdt.rm_r = + dxdt.rm_r - v1_r; 
dxdt.a     = + dxdt.a - v1_r * k.n.r * (1 + k.beta); 
 
v2_r = - S.m_r * S.e_r  *  1/10*k.b + k.u * S.rm_r; 
dxdt.m_r  = + dxdt.m_r + v2_r; 
dxdt.e_r   = + dxdt.e_r + 10 * v2_r; 
dxdt.rm_r = + dxdt.rm_r - v2_r; 
 
v3_r = + S.rm_r * k.cm * k.b_cm; 
dxdt.rm_r = + dxdt.rm_r - v3_r; 
dxdt.zm_r = + dxdt.zm_r + v3_r; 
 

 
dxdt.m_t = + dxdt.m_t + Hill(1 , k.theta.t , k.w.t , S.a); 
 
dxdt.m_t = + dxdt.m_t - k.d * S.m_t; 
 
v1_t      = + Nu(S.rm_t , k.n.t , A.gamma); 
dxdt.m_t  = + dxdt.m_t + v1_t; 
dxdt.e_r   = + dxdt.e_r  + 10 * v1_t; 
dxdt.e_t  = + dxdt.e_t + v1_t; 
dxdt.rm_t = + dxdt.rm_t - v1_t; 
dxdt.a     = + dxdt.a - v1_t * k.n.t; 
 
v2_t = - S.m_t * S.e_r  *  1/10*k.b + k.u * S.rm_t; 
dxdt.m_t  = + dxdt.m_t + v2_t; 
dxdt.e_r   = + dxdt.e_r + 10 * v2_t; 
dxdt.rm_t = + dxdt.rm_t - v2_t; 
 
v3_t = + S.rm_t * k.cm * k.b_cm; 
dxdt.rm_t = + dxdt.rm_t - v3_t; 
dxdt.zm_t = + dxdt.zm_t + v3_t; 
 

 
dxdt.m_m = + dxdt.m_m + Hill(1 , k.theta.m , k.w.m , S.a); 
 
dxdt.m_m = + dxdt.m_m - k.d * S.m_m; 
 
v1_m      = + Nu(S.rm_m , k.n.m , A.gamma); 
dxdt.m_m  = + dxdt.m_m + v1_m; 
dxdt.e_r   = + dxdt.e_r  + 10 * v1_m; 
dxdt.e_m  = + dxdt.e_m + v1_m; 
dxdt.rm_m = + dxdt.rm_m - v1_m; 
dxdt.a     = + dxdt.a - v1_m * k.n.m; 
 
v2_m = - S.m_m * S.e_r  *  1/10*k.b + k.u * S.rm_m; 
dxdt.m_m  = + dxdt.m_m + v2_m; 
dxdt.e_r   = + dxdt.e_r + 10 * v2_m; 
dxdt.rm_m = + dxdt.rm_m - v2_m; 
 
v3_m = + S.rm_m * k.cm * k.b_cm; 
dxdt.rm_m = + dxdt.rm_m - v3_m; 
dxdt.zm_m = + dxdt.zm_m + v3_m; 
 

 
dxdt.m_q = + dxdt.m_q + + Hill(1 , k.theta.q , k.w.q , S.a) * Hill(k.alpha_q , S.e_q , 1 , k.K.q); 
dxdt.m_q = + dxdt.m_q - k.d * S.m_q; 
 
v1_q      = + Nu(S.rm_q , k.n.q , A.gamma); 
dxdt.m_q  = + dxdt.m_q + v1_q; 
dxdt.e_r   = + dxdt.e_r  + 10 * v1_q; 
dxdt.e_q  = + dxdt.e_q + v1_q; 
dxdt.rm_q = + dxdt.rm_q - v1_q; 
dxdt.a     = + dxdt.a - v1_q * k.n.q; 
 
v2_q = - S.m_q * S.e_r  *  1/10*k.b + k.u * S.rm_q; 
dxdt.m_q  = + dxdt.m_q + v2_q; 
dxdt.e_r   = + dxdt.e_r + 10 * v2_q; 
dxdt.rm_q = + dxdt.rm_q - v2_q; 
 
v3_q = + S.rm_q * k.cm * k.b_cm; 
dxdt.rm_q = + dxdt.rm_q - v3_q; 
dxdt.zm_q = + dxdt.zm_q + v3_q; 
 

 
dxdt.m_p = + dxdt.m_p + Hill(1 , k.theta.p , k.w.p , S.a); 
 
dxdt.m_p = + dxdt.m_p - k.d * S.m_p; 
 
v1_p      = + Nu(S.rm_p , k.n.p , A.gamma); 
dxdt.m_p  = + dxdt.m_p + v1_p; 
dxdt.e_r   = + dxdt.e_r  + 10 * v1_p; 
dxdt.e_p  = + dxdt.e_p + v1_p; 
dxdt.rm_p = + dxdt.rm_p - v1_p; 
dxdt.a     = + dxdt.a - v1_p * k.n.p; 
 
v2_p = - S.m_p * S.e_r  *  1/10*k.b + k.u * S.rm_p; 
dxdt.m_p  = + dxdt.m_p + v2_p; 
dxdt.e_r   = + dxdt.e_r + 10 * v2_p; 
dxdt.rm_p = + dxdt.rm_p - v2_p; 
 
v3_p = + S.rm_p * k.cm * k.b_cm; 
dxdt.rm_p = + dxdt.rm_p - v3_p; 
dxdt.zm_p = + dxdt.zm_p + v3_p; 
 

end