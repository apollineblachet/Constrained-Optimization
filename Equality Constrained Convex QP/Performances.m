close all
clear
clc

% QP parameters
alpha = 0.01;
beta = 0.5;
density = 0.85;

% computing cpu time as a function of problem size
nrange = 10:10:700;
CPUTIMES = zeros(6,length(nrange));
solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace'};
for i = 1:length(nrange)
    [H,g,A,b,xtrue,lambda] =QPGenerator(nrange(i),alpha,beta,density);
    for j = 1:6
        t1 = cputime;
        [x,lambda] = EqualityQPSolver(H,g,A,b,solvers{j});
        t2 = cputime-t1;
        CPUTIMES(j,i) = t2;
    end
end

% plot results
figure
plot(nrange, CPUTIMES')
xlabel('Problem size n')
ylabel('CPU time')
legend(solvers, 'Location','northwest')