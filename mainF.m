

%Inputs:cond_numb = condition number of the matrix of the instance; qub =system syze of the matrix 2^qub; Instance; instance of the matrix used
%Output: Random matrix A and vector b: Ax=b

[A,b]=randMat_gen(cond_numb,qub,instance); 

%Inputs= A;b;qub (output and imputs from the previous function), T=Total_steps used in the methos;
%%Output: error of the solution of Ax=b

[error]=Walk_error_Herm(T,A,b,qub) 








