function [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b)
% Sparse LU
[n,m] = size(A);
[K,d] = KKTmatrix(H,g,A,b);
Ksparse = sparse(K); 
[Lsparse,Usparse,psparse] = lu(Ksparse,'vector');
z = Usparse\(Lsparse\d(psparse));
x = z(1:n);
lambda = z(n+1:length(z));
end