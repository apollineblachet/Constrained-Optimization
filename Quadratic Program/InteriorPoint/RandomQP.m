function [H,g,A,l,u,b] = RandomQP(n,alpha,density)
% RandomQP  Generates data for a random convex QP
%
%   min     0.5 x' H x + g' x
%    x
%   s.t.     A x = b
%            bl <= C x <= bu
%            l <=    x <= u
%
% Syntax: [H,g,bl,A,bu,l,u] = RandomQP(n,alpha,density)
%
%   Inputs:
%       n       number of variabels
%       alpha   regularization factor. alpha > 0
%       density density of sparse matrix. 0 < density < 1
%


m = round(n/2);

C = sprand(m,n,density);
bl = -rand(m,1) - 20;
bu = rand(m,1) + 20;

M = sprandn(n,n,density)*5;
H = M*M' + alpha*eye(n,n);
g = randn(n,1);

l = -ones(n,1)*10;
u = ones(n,1)*10;

% rewriting bounds
A = 2*m*rand(m,n); %[full(C) full(-C) eye(n,n) -eye(n,n)]'; 
b = m*m*randn(m,1);