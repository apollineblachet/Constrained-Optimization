close all
clear
clc

%% Code use for generating graphs in the report


% Allocate memory
nrange = 2:2:100;
ERRORS = zeros(1,length(nrange));
ITERATIONS = zeros(2,length(nrange));
CPUTIMES = zeros(2,length(nrange));

N = length(nrange);
i = 1;
while i <= N
    n = nrange(i);

    % Generate random LP
    [g,bl,A,bu,l,u,Aeq,beq] = RandomLP(n);

    % Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
    Abar = [full(A) full(-A) eye(n,n) -eye(n,n)]; 
    bbar = [-bl; bu; -l; u];

    % Solve LP using linprog
    t1 = cputime;
    [testsolution,fval,exit,output] = linprog(g,-Abar',bbar,Aeq',-beq);
    tquadraprog = cputime-t1;
    
    if (exit == 1) % if problem feasible

        % Solve LP using primal active-set algorithm for QPs
        t1 = cputime;
        [x,lambda,Wset,it] = QPSolverActiveSet(eye(n)*1e-3,g,Abar',-bbar,Aeq',-beq); % solve QP
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
legend('Primal Active-set solver', 'Linprog', 'Location','northwest')

figure
plot(nrange, CPUTIMES')
xlabel('Problem size n')
ylabel('CPU time')
legend('Primal Active-set solver', 'Linprog', 'Location','northwest')
