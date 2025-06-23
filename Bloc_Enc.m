function [B,dim_aux]=Bloc_Enc(M)

sz=size(M);

[U,S,V] = svd(M);

B=[U,zeros(sz(1));zeros(sz(1)),eye(sz(1))]*[S, sqrt(eye(sz(1))-S^2);sqrt(eye(sz(1))-S^2),-S]*[V',zeros(sz(1));zeros(sz(1)),eye(sz(1))];

szB=size(B);
dim_aux=szB(1)/sz(1);

end