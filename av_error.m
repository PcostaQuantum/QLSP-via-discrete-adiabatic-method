
inst=100;
q=4;

parfor k=1:inst
[M,V]=randMat_gen(10,q,k);
[error10_3q(1,k),st10_2q(1,k)]=Walk_error_Herm(34,M,V,q) %0.3466
end

parfor k=1:inst
[M,V]=randMat_gen(20,q,k);
[error20_3q(1,k),st20_2q(1,k)]=Walk_error_Herm(74,M,V,q)%0.39.8
end

parfor k=1:inst
[M,V]=randMat_gen(30,q,k);
[error30_3q(1,k),st30_2q(1,k)]=Walk_error_Herm(122,M,V,q)
end

parfor k=1:inst
[M,V]=randMat_gen(40,q,k);
[error40_3q(1,k),st40_2q(1,k)]=Walk_error_Herm(178,M,V,q)
end

parfor k=1:inst
[M,V]=randMat_gen(50,q,k);
[error50_3q(1,k),st50_2q(1,k)]=Walk_error_Herm(238,M,V,q)
end





