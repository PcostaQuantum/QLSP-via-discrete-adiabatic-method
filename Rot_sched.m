function [M]=Rot_sched(cond,p,s)

%Rs rotation that depends on the condition number of the matrix A

M=1/(sqrt((1-f(cond,p,s))^2+f(cond,p,s)^2))*[1-f(cond,p,s),f(cond,p,s);f(cond,p,s),-(1-f(cond,p,s))];

end