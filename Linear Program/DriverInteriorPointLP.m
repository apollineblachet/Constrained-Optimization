clear all 
close all
clc

%% Parameters for random QP
n=10; 

%% Generate random QP
[g,bl,A,bu,l,u,Aeq,beq]=RandomLP(n);

%% Conversion to standard form LP
[gbar, Abar, bbar] = StandardFormLP(g,bl,A,bu,l,u,Aeq,beq);

%% Linprog solution
Atest = [full(-A) full(A)]; 
btest = [-bl; bu];
xtest2 = linprog(g, Atest', btest, Aeq', -beq, l, u);

%% Interior point solution
[x,info,mu,lambda,iter] = LPippd(gbar,Abar,bbar);
x = x(n+1:2*n)-x(1:n);

%% Compare solutions
norm(x-xtest2)
