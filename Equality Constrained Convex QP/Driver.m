close all
clear
clc

% Generate random equality constrained convex QP
n = 10;
alpha = 0.01;
beta = 0.5;
density = 0.15;
[H,g,A,b,xtrue,lambda] =QPGenerator(n,alpha,beta,density);

% Solve QP with each of the 6 solvers
solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace'};
for i = 1:6
    solvers{i};
    t1 = cputime;
    [x,lambda] = EqualityQPSolver(H,g,A,b,solvers{i});
    t2 = cputime-t1;
    display([solvers{i} ' error and cpu time: '])
    err = norm(x-xtrue, 'inf')
    cpu = t2
end