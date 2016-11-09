function [ ] = Scott_plot( v , Scott)
% Make Scott figure

figure(1)

for i = 1 : 10
    
    hold all 
    plot(v.lambda(: , i) ./ log(2) , v.phi_r(: , i) , 'LineWidth' , 3);
    plot(v.lambda(: , i) ./ log(2) , v.phi_r(: , i) , 'o-');
    
end

for i = 1 : 6
    
    plot(Scott.alpha(: , i) ./ log(2) , Scott.f_r_mean(: , i) , 'o' , 'LineWidth' , 3)
        
end

ylim([0 0.6])

xlabel('Growth-rate [1/min]')
ylabel('\Phi_r')

end

