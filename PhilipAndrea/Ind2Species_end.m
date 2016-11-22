function [S] = Ind2Species_end(X) 
% gives the value to the species names from the corresponding index in the vector X 
S.m_r(: , 1) = X(: , 1); 
S.m_t(: , 1) = X(: , 2); 
S.m_m(: , 1) = X(: , 3); 
S.m_q(: , 1) = X(: , 4); 
S.m_p(: , 1) = X(: , 5); 
S.rm_r(: , 1) = X(: , 6); 
S.rm_t(: , 1) = X(: , 7); 
S.rm_m(: , 1) = X(: , 8); 
S.rm_q(: , 1) = X(: , 9); 
S.rm_p(: , 1) = X(: , 10); 
S.zm_r(: , 1) = X(: , 11); 
S.zm_t(: , 1) = X(: , 12); 
S.zm_m(: , 1) = X(: , 13); 
S.zm_q(: , 1) = X(: , 14); 
S.zm_p(: , 1) = X(: , 15); 
S.e_r(: , 1) = X(: , 16); 
S.e_t(: , 1) = X(: , 17); 
S.e_m(: , 1) = X(: , 18); 
S.e_q(: , 1) = X(: , 19); 
S.e_p(: , 1) = X(: , 20); 
S.s_i(: , 1) = X(: , 21); 
S.a(: , 1) = X(: , 22); 
end