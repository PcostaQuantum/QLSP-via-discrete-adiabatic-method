function [v]=f(k,p,s)

v=k/(k-1)*(1-(1+s*(k^(p-1)-1))^(1/(1-p)));

end

