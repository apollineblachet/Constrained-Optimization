function [g,bl,A,bu,l,u, Aeq,beq] = RandomLP(n)
% RandomQP  Generates data for a random convex LP
%
%   min      g' x
%    x
%   s.t.     Aeq' x + beq = 0
%            bl <= A' x <= bu
%            l <=    x <= u
%
% Syntax: [g,bl,A,bu,l,u, Aeq,beq]= RandomLP(n,alpha,density)
%
%   Inputs:
%       n       number of variabels
%

m = 2*n;
meq = round(0.5*n);
Aeq = randn(n,meq);
A = randn(n,m);
bl = -rand(m,1) - 20;
bu = rand(m,1) + 20;
g = randn(n,1);
l = -ones(n,1)*10;
u = ones(n,1)*10;
beq = randn(meq,1);