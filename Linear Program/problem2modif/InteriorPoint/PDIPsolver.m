function [x,lambda,iter,fval] = PDIPsolver(H,g,A,b)

n = size(H,1);

x = randn(n,1);
y = A*x - b;
k = min(y);
slack = k - 10;

z = [x; slack];
Hnew = eye(size(x,1)+1) * 0.000001;
gnew = [zeros(size(x,1),1); -1];
Anew = [A -ones(size(A,1),1)];
ynew = Anew*z - b;

iter = 0;
mu = 1;

while iter < 1000
    iter = iter + 1;
    centerparam = 0.5;
    lambda = 1./ynew;
    
    
    r1 = -Hnew*z + transpose(Anew)*lambda - gnew;
    r2 = -Anew*z + b + ynew;
    r3 = -ynew.*lambda + ones(size(ynew))*mu*centerparam;
    
    sys1 = [Hnew zeros(size(Hnew,1), size(ynew,1)) -transpose(Anew)];
    sys2 = [Anew -eye(size(ynew,1)) zeros(size(Anew,1), size(ynew,1))];
    sys3 = [zeros(size(ynew,1), size(Hnew,1)) diag(lambda) diag(ynew)];
    
    %system 16.58 fra bogen
    sys = [sys1; sys2; sys3];
    
    
    deltas = sys\[r1;r2;r3];

    dx = deltas(1:size(Hnew,1));
    dy = deltas(size(Hnew,1)+1:size(Hnew,1)+size(b));
    dlambda = deltas(size(Hnew,1)+size(b)+1:end);
   
    
    % calculating step lengths
    alpha_primal = 1;
    while min(ynew+alpha_primal*dy) < 0
        display(alpha_primal);
        alpha_primal = 0.9*alpha_primal;
    end
    
    alpha_dual = 1;
    while min(lambda+alpha_dual*dlambda) < 0
        alpha_dual = 0.9*alpha_dual;
    end
    
    %update values 
    z = z + alpha_primal * dx;
    ynew = ynew + alpha_primal * dy;
    lambda = lambda + alpha_dual * dlambda;

    mu = transpose(ynew)*lambda/size(ynew,1);
        
    
end

% initiate main problem with found feasible point.
x = z(1:n);
y = A*x - b;


%%  main optimization loop
tol = 0.0001;
iter = 0;
max_iter = 100*size(H,1);

while (mu > tol || norm(r1,'inf') > tol || norm(r2,'inf') > tol || norm(r3,'inf') > tol) && iter < max_iter
    iter = iter + 1;

    centerparam = 0.5;
    lambda = 1./y;
    
    r1 = -H*x + transpose(A)*lambda - g;
    r2 = -A*x + b + y;
    r3 = -y.*lambda + ones(size(y))*mu*centerparam;
    
    sys1 = [H zeros(size(H,1), size(y,1)) -transpose(A)];
    sys2 = [A -eye(size(y,1)) zeros(size(A,1), size(y,1))];
    sys3 = [zeros(size(y,1), size(H,1)) diag(lambda) diag(y)];
    
    sys = [sys1; sys2; sys3];
    
    
    deltas = sys\[r1;r2;r3];

    dx = deltas(1:size(H,1));
    dy = deltas(size(H,1)+1:size(H,1)+size(b));
    dlambda = deltas(size(H,1)+size(b)+1:end);
   
    
    % calculating step lengths
    alpha_primal = 1;
    while min(y+alpha_primal*dy) < 0
        %display(alpha_primal);
        alpha_primal = 0.8*alpha_primal;
    end
    
    alpha_dual = 1;
    while min(lambda+alpha_dual*dlambda) < 0
        alpha_dual = 0.8*alpha_dual;
    end
    
    %update values 
    x = x + alpha_primal * dx;
    y = y + alpha_primal * dy;
    lambda = lambda + alpha_dual * dlambda;

    mu = transpose(y)*lambda/size(y,1);
    
    fval = 0.5 *transpose(x) * H * x + transpose(g)*x;

    
end










