% Definition of rate function mu for transport and metabolism
function mu = mu(enzyme,alpha,midpoint,maxvalue,x )      
mu = enzyme * Hill(alpha,midpoint,maxvalue,x);        
end

