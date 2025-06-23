function [v]=f(k,p,s)

%k=cond(A);

v=k/(k-1)*(1-(1+s*(k^(p-1)-1))^(1/(1-p)));

end

