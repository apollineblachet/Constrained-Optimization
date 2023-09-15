function [x,lambda]=EqualityQPSolverLUdense(H,g,A,b)
% LU factorization
[n,m] = size(A);
[K,d] = KKTmatrix(H,g,A,b);
[L,U,p] = lu(K,'vector');
z = U\(L\d(p));
x = z(1:n);
lambda = z(n+1:length(z));
end