function [H,g,A,b,x,lambda] = QPGenerator(n,alpha,beta,density)
% QPGENERATOR generates a random QP together with its solution.
% The QP has the form:
%
%   min 0.5 x' H x + g' x
%    x
%   s.t.    A' x + b = 0
%     
% INPUT
%   n: size of the square matrix H
%   alpha: regularization terme to ensure that the problem is not unbounded
%   beta: term that determined the number of columns of A
%   density: density of matrices A and H
% 
% OUTPUT
%   H,g,A,b: matrices and vectors that define the QP problem
%   x: solution of the QP problem
%   lambda: Lagrandian mutipliers associated to x

m = round(beta*n);

% Generate A (must be full column rank)
r = repmat(randperm(n), 1, m); % repeat to ensure ceil(p*m*n) < length(r)
c = repmat(randperm(m), 1, n);
row_idx = r(1:ceil(density*m*n));
col_idx = c(1:ceil(density*m*n));
ran_num = rand(1, ceil(density*m*n));
A = full(sparse(row_idx, col_idx, ran_num, n, m));

% Genrate H
M = full(sprandn(n,n,density));
H = M*M' + alpha*eye(n);

% Generate x and lambda
x = randn(n,1);
lambda = randn(m,1);

% Deduce b and g using optimality conditions
b = -A'*x;
g = A*lambda - H*x;
end