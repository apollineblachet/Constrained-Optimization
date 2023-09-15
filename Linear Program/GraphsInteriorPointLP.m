close all
clear
clc

%% Code use for generating graphs in the report

% Allocate memory
nrange = 2:1:100;
ERRORS = zeros(1,length(nrange));
ITERATIONS = zeros(2,length(nrange));
CPUTIMES = zeros(2,length(nrange));

N = length(nrange);
i = 1;
while i <= N
    n = nrange(i);

    % Generate random LP
    [g,bl,A,bu,l,u,Aeq,beq]=RandomLP(n);

    % Conversion to standard form LP
    [gbar, Abar, bbar] = StandardFormLP(g,bl,A,bu,l,u,Aeq,beq);

    % Solve LP using linprog
    Atest = [full(-A) full(A)]; 
    btest = [-bl; bu];
    t1 = cputime;
    [testsolution,fval,exit,output] = linprog(g, Atest', btest, Aeq', -beq, l, u);
    tlinprog = cputime-t1;
    
    if (exit == 1) % if problem feasible

        % Solve QP using primal active-set algorithm
        t1 = cputime;
        [x,info,mu,lambda,iter] = LPippd(gbar,Abar,bbar);
        x = x(n+1:2*n)-x(1:n);
        t2 = cputime-t1;
        
        % Store results
        ERRORS(1,i) = norm(x-testsolution,'inf');
        ITERATIONS(1:2,i) = [iter, output.iterations];
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
legend('Primal-dual Interior point', 'Linprog', 'Location','northwest')

figure
plot(nrange, CPUTIMES')
xlabel('Problem size n')
ylabel('CPU time')
legend('Primal-dual Interior point', 'Linprog', 'Location','northwest')
