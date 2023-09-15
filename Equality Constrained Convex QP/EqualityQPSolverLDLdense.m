function [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b)
% LDL factorization
[n,m] = size(A);
[K,d] = KKTmatrix(H,g,A,b);
[Ld,Dd,pd] = ldl(K,'vector');
z = zeros(n+m,1);
z(pd) = Ld' \ (Dd \ ( Ld \ d(pd)));
x = z(1:n);
lambda = z(n+1:length(z));
end