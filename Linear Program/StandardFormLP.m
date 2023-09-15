function [gbar, Abar, bbar] = StandardFormLP(g,bl,A,bu,l,u,Aeq,beq)
%STANDARDFORMLP Convert LP from assignement form to standard form
%   See Linear Programming, Lecture 07A - The Revised Simplex Algorithm,
%   slide 11
%
% FROM:
%   min     g' x
%    x
%   s.t.     Aeq'x + beq = 0
%            bl <= A' x <= bu
%            l <=    x <= u
%
% TO:
%   min     gbar' x
%    x
%   s.t.     Abar'x + bbar = 0
%            x >= 0


n = size(A,1);
m = size(A,2);
neq = size(Aeq,1);
meq = size(Aeq,2);
gbar = [-g; g; zeros(m,1); zeros(m,1); zeros(n,1); zeros(n,1)];
bbar = [-beq; bl ; bu; l; u];
Abar = [-Aeq',    Aeq',    zeros(meq,m), zeros(meq,m), zeros(meq,n), zeros(meq,n);
        -A' ,     A',      -eye(m),      zeros(m,m),   zeros(m,n),   zeros(m,n); 
        -A' ,     A' ,     zeros(m,m),   eye(m),       zeros(m,n),   zeros(m,n); 
        -eye(n),  eye(n),  zeros(n,m),   zeros(n,m),   -eye(n),      zeros(n,n); 
        -eye(n),  eye(n),  zeros(n,m),   zeros(n,m),   zeros(n,n),   eye(n)];

end

