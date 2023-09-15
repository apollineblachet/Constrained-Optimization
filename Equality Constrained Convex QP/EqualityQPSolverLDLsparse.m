function [x,lambda]=EqualityQPSolverLDLsparse(H,g,A,b)
% Sparse LDL
[n,m] = size(A);
[K,d] = KKTmatrix(H,g,A,b);
Ksparse = sparse(K); 
[Ldsparse,Ddsparse,pdsparse] = ldl(Ksparse,'vector');
z = zeros(n+m,1);
z(pdsparse) = Ldsparse' \ (Ddsparse \ ( Ldsparse \ d(pdsparse)));
x = z(1:n);
lambda = z(n+1:length(z));
end