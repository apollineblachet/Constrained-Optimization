clear all
close all
clc

%% Parameters for random LP
n=30; 

%% Generate random LP
[g,bl,A,bu,l,u,Aeq,beq] = RandomLP(n);

%% Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
Abar = [full(A) full(-A) eye(n,n) -eye(n,n)]; 
bbar = [-bl; bu; -l; u];

%% Solve LP using linprog
testsolution = linprog(g,-Abar',bbar);

%% solve LP using Primal Dual Interior Point method for QPs
[PDIPx,PDIPl,iter,fval] = QPSolverPrimalDualIP(eye(n)*1e-3,g,Abar',-bbar,[],[]);

%% Comparison
norm(testsolution - PDIPx,'inf')



 
