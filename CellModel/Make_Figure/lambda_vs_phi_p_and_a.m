

c  = 250;
va = 5;
Ka = 5;
k2 = 1;

Rate_function = @(x,y) (va * c * x) ./ (Ka ./ y + 1 + va) * 1 / k2; 

Ia_i = Rate_function(v.phi_t , v.a);

%%
figure(1)

hold all

plot( 60 ./ (log(2) ./ v.lambda) , v.phi_t * 150 , 'Color' , 'blue' , 'Marker' , 'd' , 'LineWidth' , 3)
plot( 60 ./ (log(2) ./ v.lambda) , v.a           , 'Color' , 'red'  , 'Marker' , 'd' , 'LineWidth' , 3)
plot( 60 ./ (log(2) ./ v.lambda) , Ia_i          , 'Color' , 'black', 'Marker' , 'd' , 'LineWidth' , 3)

legend('P protein' , 'a' , '[I^+]' , 'Location' , 'NorthWest')

x_label = xlabel('Growth-Rate \lambda (1/h)');
y_label = ylabel('Concentration [a.u.]');

set(y_label ,'FontSize',16);
set(x_label ,'FontSize',16);

print(gcf , 'Figures/lambda_P_and_a_2.eps' , '-dpsc2')