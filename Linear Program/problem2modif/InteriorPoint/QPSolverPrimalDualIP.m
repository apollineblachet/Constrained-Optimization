function [x,lambda,iter,fval] = QPSolverPrimalDualIP(H,g,A,b,Aeq,beq)

sizeA = size(A,1);
sizeAeq = size(Aeq,1);


if sizeA < 1 
    
    [n,m] = size(A);
    K = [H, A; -A', zeros(m,m)];
    
    % making sure the KKT matrix is dense:
    K = full(K);
    
    d = [-g;-b];
    
    z = zeros(n+m,1);
    [L,D,p] = ldl(K,'lower','vector');
    z(p) = L' \ ( D \ ( L \ d(p) ) );
    
    x = z(1:size(H,1));
    lambda = z(size(H,1)+1:end);

    iter = 0;

elseif sizeAeq < 1
    [x,lambda,iter,fval] = PDIPsolver(H,g,A,b);
else
    [x,lambda,iter,fval] = PDIPsolver_wEQ(H,g,A,b,Aeq,beq);
end