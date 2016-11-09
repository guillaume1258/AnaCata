 function y = Hill(alpha,midpoint,maxvalue,x)
 % Hill function
 
            y = ...
                maxvalue * (1/ (1 + (midpoint/x)^alpha));
end
