close all
clear
clc

%% Code use for generating graphs in the report

% Allocate memory
nrange = 2:1:100;
ERRORS = zeros(1,length(nrange));
ITERATIONS = zeros(3,length(nrange));
CPUTIMES = zeros(2,length(nrange));

N = length(nrange);
i = 1;
while i <= N
    n = nrange(i);

    % Generate random LP
    [g,bl,A,bu,l,u,Aeq,beq]=RandomLP(n);

    % Conversion to standard form LP
    [g, A, b] = StandardFormLP(g,bl,A,bu,l,u,Aeq,beq);

    % Solve LP using linprog
    t1 = cputime;
    [testsolution,fval,exit,output] = linprog(g, [], [], A, b, zeros(size(g,1),1), []);
    tlinprog = cputime-t1;
    
    if (exit == 1) % if problem feasible

        % Solve QP using primal active-set algorithm
        t1 = cputime;
        [x, f_opt, n1, n2] = RevisedSimplexComplete(A,b,g);
        t2 = cputime-t1;
        
        % Store results
        ERRORS(1,i) = norm(x-testsolution,'inf');
        ITERATIONS(1:3,i) = [n1, n2, output.iterations];
        CPUTIMES(1:2,i) = [t2,tlinprog];

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
legend('Simplex Phase I', 'Simplex Phase II', 'Linprog', 'Location','northwest')

figure
plot(nrange, CPUTIMES')
xlabel('Problem size n')
ylabel('CPU time')
legend('Simplex', 'Linprog', 'Location','northwest')
