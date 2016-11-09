function [ ] = Scott( v )
% Make Scott figure
figure(1)
subplot(2 , 2 , 1)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.phi_r(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.phi_r(: , i) , 'o-');
end

ylim([0 0.6])

xlabel('\lambda')
ylabel('\Phi_r')

subplot(2 , 2 , 2)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.phi_t(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.phi_t(: , i) , 'o-');
end

xlabel('\lambda')
ylabel('\Phi_t (= \Phi_m)')

subplot(2 , 2 , 3)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.phi_q(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.phi_q(: , i) , 'o-');
end

xlabel('\lambda')
ylabel('\Phi_q')

subplot(2 , 2 , 4)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.a(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.a(: , i) , 'o-');
end

xlabel('\lambda')
ylabel('a')


saveas(gcf , 'Figures/physiology_cm.pdf' , 'pdf')

%% plot total mRNA content as well as the processivity per ribosome.
close all

figure(2)
subplot(2 , 2 , 1)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.mRNA_tot(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.mRNA_tot(: , i) , 'o-');
end

xlabel('\lambda')
ylabel('Total mRNA molecules')

subplot(2 , 2 , 2)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.zm_rm_ratio(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.zm_rm_ratio(: , i) , 'o-');
end

xlabel('\lambda')
ylabel('% of zombie complexes')

subplot(2 , 2 , 3)
for i = 1 : 6
    hold all 
    plot(v.lambda(: , i) , v.gamma(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) , v.gamma(: , i) , 'o-');
end

xlabel('\lambda')
ylabel('\gamma')

saveas(gcf , 'Figures/physiology_cm_2.pdf' , 'pdf')

end

