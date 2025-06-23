function [Aproxs,steps]=Walk_error_Herm(T,A,V,n)


Bl_A=Bloc_Enc(A);
Id=eye(2);
ket0=[1;0];
ket1=[0;1];


Proj0=kron(ket0,ket0');
Proj1=kron(ket1,ket1');
Had=1/sqrt(2)*[1,1;1,-1];
X=[0,1;1,0];
Z=[1,0;0,-1];
borb01=kron(ket0,ket1');
borb10=kron(ket1,ket0');

Had_sys = kron(kron(kron(kron(kron(Id,Id),Had),Id),Id),eye(2^n));
X_syst =  kron(kron(kron(kron(kron(Id,Id),Id),X),Id),eye(2^n));
Ub_syst = kron(kron(kron(kron(kron(Id,Id),Id),Id),Id),UPrep_b(V,n));



Cont_UA = kron(kron(kron(kron(kron(Z,Proj0),Id),Id),Id),eye(2^n)) + kron(kron(kron(kron(borb01,Proj1),Id),Id),Bl_A)+...
kron(kron(kron(kron(borb10,Proj1),Id),Id),Bl_A');

CUQb0= kron(kron(kron(kron(kron(Id,Id),Id),Proj1),Id),eye(2^n))+Ub_syst*(kron(kron(kron(kron(kron(Id,Id),Proj0),Proj0),Id),eye(2^n))+...
    kron(kron(kron(kron(kron(Id,Id),Proj1),Proj0),Id),eye(2^n))-kron(kron(kron(kron(kron(Proj0,Id),Proj1),Proj0),Id),2*binquant(0,n)*binquant(0,n)'))*Ub_syst';

CUQb1= kron(kron(kron(kron(kron(Id,Id),Id),Proj0),Id),eye(2^n))+Ub_syst*(kron(kron(kron(kron(kron(Id,Id),Proj0),Proj1),Id),eye(2^n))+...
  kron(kron(kron(kron(kron(Id,Id),Proj1),Proj1),Id),eye(2^n))-kron(kron(kron(kron(kron(Proj0,Id),Proj1),Proj1),Id),2*binquant(0,n)*binquant(0,n)'))*Ub_syst';

C_Had=kron(kron(kron(kron(kron(Id,Id),Id),Proj0),Id),eye(2^n)) + kron(kron(kron(kron(kron(Id,Had),Id),Proj1),Id),eye(2^n));


Prj=kron(kron(kron(kron(kron(Id,ket0),ket0),Id),ket0),eye(2^n))';

Ref= 2*kron(kron(kron(kron(kron(Id,Proj0),Proj0),Id),Proj0),eye(2^n)) - eye(2^(n+5));

b=UPrep_b(V,n)*binquant(0,n);
x_exact = A\b;  
x_exact_n = x_exact/norm(x_exact);
x_exact_n2=kron(kron([0;1],[1;0]),x_exact_n);

In=kron(kron(kron(kron(kron(ket0,ket0),ket0),ket0),ket0),UPrep_b(V,n)*binquant(0,n));
for k=1:T
   Rot_Sc=Rot_sched(cond(A),1.4,(k)/j);
   C_sched=kron(kron(kron(kron(kron(Id,Rot_Sc),Id),Proj0),Id),eye(2^n)) + kron(kron(kron(kron(kron(Id,Id),Id),Proj1),Id),eye(2^n));
   In=Had_sys*In;
   In=CUQb1*In;
   In=C_sched*In;
   In=C_Had*In;
   In=Cont_UA*In;
   In=X_syst*In;
   In=C_Had*In;
   In=C_sched*In;
   In=CUQb1*In;
   In=Had_sys*In;
   In = Ref*In;
  end
  w=Prj*In;

Aproxs=norm(x_exact_n2-w/norm(w));
steps=norm(w);

end



