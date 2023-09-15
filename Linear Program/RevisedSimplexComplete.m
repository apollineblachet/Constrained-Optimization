function [x_opt, f_opt, n1, n2] = RevisedSimplexComplete(A,b,g)

%% Implementation of the revised Simplex using Phase I method for finfing and initial point. Solves a linear
% programming problem of the form
%
%   min g'*x
%   s.t. Ax  = b
%         x >= 0
%
% The function input parameters are the following:
%     A: The constraint matrix 
%     b: The rhs vector 
%     g: The vector of cost coefficients 
% 
% The function returns:
%     x_opt: Decision variable values at the optimal solution  
%     f_opt: Objective function value at the optimal solution
%     n1: number of iterations of initialization
%     n2: number of iterations in main algorithm


%% Phase I (see Lecture 07A - The Revised Simplex Algorithm, slide 24)

% Construct phase I problem in standard form
m = size(A,1);
n = size(A,2);
gPI = [zeros(n,1); 1; zeros(m,1); zeros(m,1)];
API = [A, ones(m,1), -eye(m), zeros(m,m);
    -A, ones(m,1), zeros(m,m), -eye(m)];
bPI = [b;-b];

% Initial feasible point and associated basic set (its size should be 2m)
t = max(abs(b));
x0 = [zeros(n,1); max(abs(b)); t*ones(m,1)-b; t*ones(m,1)+b];
N = find(abs(x0) < 1e-9);
N = N(1:n+1);
B = 1:2*m+n+1;
B(N) = [];

% Solve
[x0_opt, f0_opt, n1] = RevisedSimplex(API,bPI,gPI,B,x0);

%% Phase II (see Lecture 07A - The Revised Simplex Algorithm, slide 26)
x0 = x0_opt(1:n); % Report the solution in terms of the original problem.
if (norm(f0_opt)<1e-9)
    N = find(abs(x0) < 1e-9);
    N = N(1:n-m);
    B = 1:n;
    B(N) = [];
    [x_opt, f_opt, n2] = RevisedSimplex(A,b,g,B,x0);
else
    display('Problem infeasible')
    f_opt = -inf;  % Optimal objective function cost = -inf
    x_opt = [];    % Produce empty vector []
    return         % Break from the function
end

end