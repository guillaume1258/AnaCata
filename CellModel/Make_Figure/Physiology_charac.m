figure(1)
subplot(2 , 2 , 1)

pcolor(v.n_s , v.cm , v.phi_r)
%set(gca , 'XTick' , log(2) ./ [v.lambda])
%
subplot(2 , 2 , 2)

pcolor(v.n_s , v.cm , v.phi_t)




%
subplot(2 , 2 , 3)
pcolor(v.n_s , v.cm , v.phi_q)




%
subplot(2 , 2 , 4)
pcolor(v.n_s , v.cm , v.a)

figure(2)
pcolor(v.n_s , v.cm , log(log(2) ./ v.lambda))