close all
clear
clc

%% Code use for generating graphs in the report


% Allocate memory
nrange = 2:2:50;
ERRORS = zeros(1,length(nrange));
ITERATIONS = zeros(2,length(nrange));
CPUTIMES = zeros(2,length(nrange));

N = length(nrange);
i = 1;
while i <= N
    n = nrange(i);

    % Generate random QP
    [g,bl,A,bu,l,u,Aeq,beq] = RandomLP(n);

    %% Convert bl <= A' x <= bu, l <= x <= u to c(x) = Abar' x + bbar >= 0
    Abar = [full(A) full(-A) eye(n,n) -eye(n,n)]; 
    bbar = [-bl; bu; -l; u];

    % Solve QP using quadraprog
    t1 = cputime;
    [testx,fval,exitflag,output] = linprog(g,-Abar',bbar);
    tquadraprog = cputime-t1;
    
    if (exitflag == 1) % if problem feasible

        % Solve QP using primal active-set algorithm
        t1 = cputime;
        [PDIPx,PDIPl,iter,fval] = QPSolverPrimalDualIP(eye(n)*1e-3,g,Abar',-bbar,[],[]); % solve QP
        t2 = cputime-t1;
        
        % Store results
        ERRORS(1,i) = norm(testx - PDIPx,'inf');
        ITERATIONS(1:2,i) = [iter, output.iterations];
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
legend('Primal-dual Interior point', 'Linprog', 'Location','northwest')

figure
plot(nrange, CPUTIMES')
xlabel('Problem size n')
ylabel('CPU time')
legend('Primal-dual Interior point', 'Linprog', 'Location','northwest')
