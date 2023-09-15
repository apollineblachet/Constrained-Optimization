function [x,info,mu,lambda,iter] = LPippd(g,A,b)

% LPIPPD   Primal-Dual Interior-Point LP Solver
%
%          min  g'*x
%           x
%          s.t. A x  = b      (Lagrange multiplier: mu)
%                 x >= 0      (Lagrange multiplier: lamba)
%
% Syntax: [x,info,mu,lambda,iter] = LPippd(g,A,b)
%
%         info = true   : Converged
%              = false  : Not Converged

%%
[m,n]=size(A);

maxit = 100;
tolL = 1.0e-9;
tolA = 1.0e-9;
tols = 1.0e-9;

eta = 0.99;

%% STARTING POINT p.410 Nocedal

x = A'*inv(A*A')*b;
mu = (A*A')\A*g;
lambda = g - A'*mu;

ddx = max(-3/2*min(x),0);
ddlambda = max(-3/2*min(lambda),0);
x = x +ddx*ones(n,1);
lambda = lambda +ddlambda*ones(n,1);

ddx = x'*lambda/sum(lambda)/2;
ddlambda = x'*lambda/sum(x)/2;
x = x +ddx*ones(n,1);
lambda = lambda +ddlambda*ones(n,1);

%%

% Compute residuals
rL = g - A'*mu - lambda;    % Lagrangian gradient
rA = A*x - b;               % Equality Constraint
rC = x.*lambda;             % Complementarity
s = sum(rC)/n;              % Duality gap

% Converged
Converged = (abs(s) <= tols);
%%        
iter = 0;
while ~Converged && (iter<maxit)
    iter = iter+1;
    
    % ====================================================================
    % Form and Factorize Hessian Matrix
    % ====================================================================
    xdivlambda = x./lambda;
    H = A*diag(xdivlambda)*A'; % eq 62 slide 21 week 7
    L = chol(H,'lower');
   
    % ====================================================================
    % Affine Step
    % ====================================================================
    % Solve
    tmp = (x.*rL + rC)./lambda; % eq 62 slide 21 week 7
    rhs = -rA + A*tmp; % eq 62 slide 21 week 7 
    
    dmu = L'\(L\rhs); % eq 62 slide 21 week 7
    dx = xdivlambda.*(A'*dmu) - tmp; % eq 60 slide 21 week 7
    dlambda = -(rC+lambda.*dx)./x; % eq 55 slide 20 week 7
    
    % Step length
    idx = find(dx < 0.0);
    alpha = min([1.0; -x(idx,1)./dx(idx,1)]);
    
    idx = find(dlambda < 0.0);
    beta = min([1.0; -lambda(idx,1)./dlambda(idx,1)]);
    
    % ====================================================================
    % Center Parameter
    % ====================================================================
    xAff = x + alpha*dx;
    lambdaAff = lambda + beta*dlambda;
    sAff = sum(xAff.*lambdaAff)/n;
    
    sigma = (sAff/s)^3;
    tau = sigma*s;    

    % ====================================================================
    % Center-Corrector Step
    % ====================================================================
    rC = rC + dx.*dlambda - tau; % change rc and then do the same as before
    
    tmp = (x.*rL + rC)./lambda;
    rhs = -rA + A*tmp;
    
    dmu = L'\(L\rhs);
    dx = xdivlambda.*(A'*dmu) - tmp;
    dlambda = -(rC+lambda.*dx)./x;
    
    % Step length
    idx = find(dx < 0.0);
    alpha = min([1.0; -x(idx,1)./dx(idx,1)]);
    
    idx = find(dlambda < 0.0);
    beta = min([1.0; -lambda(idx,1)./dlambda(idx,1)]);

    % ====================================================================
    % Take step 
    % ====================================================================
    x = x + (eta*alpha)*dx;
    mu = mu + (eta*beta)*dmu;
    lambda = lambda + (eta*beta)*dlambda;

 
    
    % ====================================================================
    % Residuals and Convergence
    % ====================================================================
    % Compute residuals
    rL = g - A'*mu - lambda;    % Lagrangian gradient
    rA = A*x - b;               % Equality Constraint
    rC = x.*lambda;             % Complementarity
    s = sum(rC)/n;              % Duality gap

    % Converged (just check duality gap otherwise rA increase
    % again while waiting rL or s to decrease => convergence is never
    % reached)
    Converged =  (abs(s) <= tols);
end

%%
% Return solution
info = Converged;
if ~Converged
    x=[];
    mu=[];
    lambda=[];
end
    