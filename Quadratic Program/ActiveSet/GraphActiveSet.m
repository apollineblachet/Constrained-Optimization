close all
clear
clc

%% Code use for generating graphs in the report

% QP parameters
alpha = 0.01;
density = 0.15;

% computing cpu time as a function of problem size
nrange = 2:100;
ERRORS = zeros(1,length(nrange));
ITERATIONS = zeros(2,length(nrange));
CPUTIMES = zeros(2,length(nrange));

N = length(nrange);
i = 1;
while i <= N
    n = nrange(i);

    % Generate random QP
    [H,g,bl,A,bu,l,u,Aeq,beq] = RandomQP(n,alpha,density);

    % Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
    Abar = [full(A) full(-A) eye(n,n) -eye(n,n)]; 
    bbar = [-bl; bu; -l; u];

    % Solve QP using quadraprog
    t1 = cputime;
    [testsolution,fval,exit,output] = quadprog(H,g,-Abar',bbar,Aeq',-beq);
    tquadraprog = cputime-t1;
    
    if (exit == 1) % if problem feasible

        % Solve QP using primal active-set algorithm
        t1 = cputime;
        [x,lambda,Wset,it] = QPSolverActiveSet(H,g,Abar',-bbar,Aeq',-beq); % solve QP
        t2 = cputime-t1;
        
        % Store results
        ERRORS(1,i) = norm(x-testsolution,'inf');
        ITERATIONS(1:2,i) = [it, output.iterations];
        CPUTIMES(1:2,i) = [t2,tquadraprog];

        i = i+1;
    end
end

% plot results
figure
plot(nrange, ERRORS')
xlabel('Problem size n')
ylabel('Norm inf of the error')

figure
plot(nrange, ITERATIONS')
xlabel('Problem size n')
ylabel('Number of itérations')
legend('Primal Active-set solver', 'Quadraprog', 'Location','northwest')

figure
plot(nrange, CPUTIMES')
xlabel('Problem size n')
ylabel('CPU time')
legend('Primal Active-set solver', 'Quadraprog', 'Location','northwest')
