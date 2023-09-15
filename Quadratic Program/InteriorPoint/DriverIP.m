clear all
close all
clc

%% Parameters for random QP
n=30; 
alpha=0.01; 
density=0.15; 

%% Generate random QP
[H,g,A,l,u,b] = RandomQP(n,alpha,density);

%% calulate solution using quadprog
[testx,fval,exitflag,output,testlambda] = quadprog(H,g,-A,-b);

%% solve using Primal Dual Interior Point method
[PDIPx,PDIPl,iter,fval] = QPSolverPrimalDualIP(H,g,A,b,[],[]);

%% Comparison
norm(testx - PDIPx,'inf')



 
