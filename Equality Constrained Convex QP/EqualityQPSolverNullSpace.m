function [x,lambda]=EqualityQPSolverNullSpace(H,g,A,b)
% Null space factorization
[n,m] = size(A);
[Qnull,Rbarnull] = qr(A);
m1null = size(Rbarnull,2);
Q1null = Qnull(:,1:m1null);
Q2null = Qnull(:,m1null+1:n);
Rnull = Rbarnull(1:m1null,1:m1null);

xYnull = Rnull' \ (-b); % minus b because of convention
HQnull = Q2null'*H*Q2null;
bZnull = -Q2null'*(H*Q1null*xYnull + g);
xZnull = HQnull \ bZnull;
x = Q1null*xYnull + Q2null*xZnull;
bYnull = Q1null'*(H*x+g);
lambda = Rnull \ bYnull;
end