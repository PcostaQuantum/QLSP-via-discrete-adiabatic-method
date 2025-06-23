function [ket]=binquant(v,n)


b=int2bit(v,n);

ket0=[1;0];
ket1=[0;1];

if b(n)==1
    In=ket1;
else
    In=ket0;
end

for j=1:n-1
    if b(n-j)==1
        In=kron(In,ket1);
    else
        In=kron(In,ket0);
    end
end

ket=In;
end
