clear all
close all
clc

%% Parameters for random QP
n=10; 
alpha=0.01; 
density=0.15; 

%% Generate random QP
[H,g,bl,A,bu,l,u,Aeq,beq] = RandomQP(n,alpha,density); 

%% Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
Abar = [full(A) full(-A) eye(n,n) -eye(n,n)]; 
bbar = [-bl; bu; -l; u];

%% Solve QP using quadraprog
testsolution = quadprog(H,g,-Abar',bbar,Aeq',-beq);

%% Solve QP using active set algorithm
[x,lambda,Wset,it] = QPSolverActiveSet(H,g,Abar',-bbar,Aeq',-beq); % solve QP
    
%% Compare solutions
norm(x-testsolution,'inf')
