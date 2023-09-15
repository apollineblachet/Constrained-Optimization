function [x,lambda]=EqualityQPSolverRangeSpace(H,g,A,b)
% Range space factorization
[Lrange] = chol(H,'lower');
vrange = Lrange' \ ( Lrange\g);
Hrange_a = A'* (Lrange' \ ( Lrange\A));
[Lrange_a] = chol(Hrange_a,'lower');
brange_1 = -b + A'*vrange;  % minus b because of convention
lambda = Lrange_a' \ ( Lrange_a\brange_1 );
brange_2 = A*lambda - g;
x = Lrange' \ ( Lrange\brange_2 );
end