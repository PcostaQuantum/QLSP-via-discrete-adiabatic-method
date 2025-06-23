function [Ub]= UPrep_b(b)

sz=size(b);
v(1,1)=sqrt((b(1,1)+1)/2);
for k=2:sz(1)
    v(k,1)=b(k,1)/(2*v(1,1));
end

Ub=2*v*v'-eye(sz(1));

end
