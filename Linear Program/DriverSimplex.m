clear all 
close all
clc

%% Parameters for random QP
n=10; 

%% Generate random QP
[g,bl,A,bu,l,u,Aeq,beq]=RandomLP(n);

%% Conversion to standard form LP
[g, A, b] = StandardFormLP(g,bl,A,bu,l,u,Aeq,beq);

%% Linprog solution
solstform = linprog(g, [], [], A, b, zeros(size(g,1),1), []);

%% Simplex solution
[x_opt, f_opt, n1, n2] = RevisedSimplexComplete(A,b,g);

%% Compare solutions
norm(solstform-x_opt)
