function [x,lambda,iter,fval] = PDIPsolver_wEQ(H,g,A,b,Aeq,beq)

n = size(H,1);

x = randn(n,1);


s = A*x - b;
k = min(s);
slack = k - 10;

z = [x; slack];
Hnew = eye(size(x,1)+1) * 0.000001;
gnew = [zeros(size(x,1),1); -1];
Aeqnew = [Aeq -ones(size(Aeq,1),1)];
Anew = [A -ones(size(A,1),1)];
ynew = Aeqnew*z - beq;
snew = Anew*z - b;
lambda = ones(size(snew));


mu = 1;

iter = 0;
max_iter = 1000*size(H,1);

while iter < max_iter
    iter = iter + 1;

    centerparam = 0.5;

    mu = transpose(snew)*lambda/size(snew,1);

    r1 = -Hnew*z + Aeqnew'*ynew - gnew + Anew'*lambda ;
    r2 = Aeqnew*z - beq;
    r3 = Anew*z - b - snew;  
    r4 = -snew.*lambda + ones(size(snew))*mu*centerparam;
    
    sys1 = [Hnew -Aeqnew' -Anew' zeros(size(Hnew,1),size(snew,2))];
    sys2 = [-Aeqnew zeros(size(Aeqnew,1), size(Aeqnew,1)+size(Anew,1)+size(snew,1))];
    sys3 = [-Anew zeros(size(Anew,1),size(Aeqnew,1) + size(Anew,1)), eye(size(snew,1))];
    sys4 = [zeros(size(snew,2), size(Hnew,2)+ size(Aeqnew,1)) diag(lambda) diag(snew)];
    
    sys = [sys1; sys2; sys3; sys4];
    
    deltas = sys\[r1;r2;r3;r4];

    dz = deltas(1:size(Hnew,1));
    dy = deltas(size(Hnew,1)+1:size(Hnew,1)+size(beq));
    dlambda = deltas(size(Hnew,1)+size(beq)+1:size(Hnew,1)+size(beq)+size(b));
    ds = deltas(size(Hnew,1)+size(beq)+1+size(b):end);
   
   
    if ds < 0 || dlambda < 0 || snew < 0
        alpha_primal = 1;
        alpha_dual = 1;
    else
        % calculating step lengths
        alpha_primal = 1;
        while min(s+alpha_primal*ds) < 0
            alpha_primal = 0.90*alpha_primal;
         
        end
        
        alpha_dual = 1;
        while min(lambda+alpha_dual*dlambda) < 0
            alpha_dual = 0.90*alpha_dual;
        end
    end
    

    %update values 
    z = z + alpha_primal * dz;
    ynew = ynew + alpha_primal * dy;
    snew = snew + alpha_primal * ds;
    lambda = lambda + alpha_dual * dlambda;

    mu = transpose(snew)*lambda/size(snew,1);

end


% initiate main problem with found feasible point.
x =  z(1:n);
y = Aeq*x - beq;
s = A*x - b;
lambda = ones(size(s));


mu = 1;

%%  main optimization loop

iter = 0;
max_iter = 10000;
tol = 0.00001;
r1 = 1;
r2 = 1;
r3 = 1;

while (mu > tol || norm(r1,'inf') > tol || norm(r2,'inf') > tol || norm(r3,'inf') > tol) && iter < max_iter
    
    iter = iter + 1;
    centerparam = 0.5;
    

    r1 = -H*x + Aeq'*y - g + A'*lambda ;
    r2 = Aeq*x - beq;
    r3 = A*x - b - s;  
    r4 = -s.*lambda + ones(size(s))*mu*centerparam;
    
    sys1 = [H -Aeq' -A' zeros(size(H,1),size(s,2))];
    sys2 = [-Aeq zeros(size(Aeq,1), size(Aeq,1)+size(A,1)+size(s,1))];
    sys3 = [-A zeros(size(A,1),size(Aeq,1) + size(A,1)), eye(size(s,1))];
    sys4 = [zeros(size(s,2), size(H,2)+ size(Aeq,1)) diag(lambda) diag(s)];
    
    sys = [sys1; sys2; sys3; sys4];
    
    deltas = sys\[r1;r2;r3;r4];

    dx = deltas(1:size(H,1));
    dy = deltas(size(H,1)+1:size(H,1)+size(beq));
    dlambda = deltas(size(H,1)+size(beq)+1:size(H,1)+size(beq)+size(b));
    ds = deltas(size(H,1)+size(beq)+1+size(b):end);


    if ds < 0 || dlambda < 0 || snew < 0
        display("Stuck");
    end
    % calculating step lengths
    alpha_primal = 1;
    while min(s+alpha_primal*ds) < 0
        alpha_primal = 0.90*alpha_primal;
      
    end
    
    alpha_dual = 1;
    while min(lambda+alpha_dual*dlambda) < 0
        alpha_dual = 0.90*alpha_dual;
     
    end
    
    mu = transpose(s)*lambda/size(s,1);

    %update values 
    x = x + alpha_primal * dx;
    y = y + alpha_primal * dy;
    s = s + alpha_primal * ds;
    lambda = lambda + alpha_dual * dlambda;


    muaff = transpose(s)*lambda/size(s,1);

    centerparam = (muaff/mu)^3;

    fval = 0.5 *transpose(x) * H * x + transpose(g)*x;
end




















