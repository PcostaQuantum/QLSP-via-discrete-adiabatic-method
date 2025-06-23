function [M]=Rot_sched(cond,p,s)

M=1/(sqrt((1-f(cond,p,s))^2+f(cond,p,s)^2))*[1-f(cond,p,s),f(cond,p,s);f(cond,p,s),-(1-f(cond,p,s))];

end
