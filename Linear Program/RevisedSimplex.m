function [x_opt, f_opt, n] = RevisedSimplex(A,b,c,C,x)

%% Implementation of the revised Simplex. Require initial point and initial basic set. Solves a linear
% programming problem of the form
%
%   min c'*x
%   s.t. Ax  = b
%         x >= 0
%
% The function input parameters are the following:
%     A: The constraint matrix 
%     b: The rhs vector 
%     c: The vector of cost coefficients 
%     C: The indices of the basic variables corresponding to an
%        initial basic feasible solution
%     x0: An initial feasible point
% 
% The function returns:
%     x_opt: Decision variable values at the optimal solution  
%     f_opt: Objective function value at the optimal solution
%     n: number of iterations 
%


    %% Initialization phase 
    % Create the initial Basis matrix, compute its inverse and    
    % compute the inital basic feasible solution
    [mA,nA]=size(A);
    B=A(:,C);
    invB = inv(B);
    x(C) = invB*b; % (13.18)


    %% Iteration phase
    n_max = 10000;      % At most n_max iterations
    for n = 1:n_max    % Main loop

        % Compute the vector of reduced costs c_r 
        c_B = c(C);         % Basic variable costs
        p = (c_B'*invB)';   % Dual variables (13.20)
        c_r = c' - p'*A;    % Vector of reduced costs (13.21)

        % Check if the solution is optimal. If optimal, use 

        [M,J] = min(c_r); % Find indices with negative reduced costs

        if (M > -1e-9) % optimal solution found
            f_opt = c'*x;
            x_opt = x;
            return;
        end

        % Choose the entering variable
        j_in = J;

        % Compute the vector u (i.e., the reverse of the basic directions) 

        u = invB*A(:,j_in);
        I = find(u > 1e-9);

        if (isempty(I)) % unbounded problem, no solution
           f_opt = -inf;  % Optimal objective function cost = -inf
           x_opt = [];    % Produce empty vector []
           return         % Break from the function
        end

        % Compute the optimal step length theta

        theta = min(x(C(I))./u(I));

        L = find(x(C)./u == theta); % Find all indices with ratio theta

        % Select the exiting variable

        l_out = L(1);

        % Move to the adjacent solution 

        x(C) = x(C) - theta*u;

        % Value of the entering variable is theta
        x(j_in) = theta;


        % Update the set of basic indices C 

        C(l_out) = j_in;

        % Compute the new inverse basis B^-1 by performing elementary row 
        % operations on [B^-1 u] (pivot row index is l_out). The vector u is trans-
        % formed into a unit vector with u(l_out) = 1 and u(i) = 0 for
        % other i.

        M=horzcat(u, invB);
        [f g]=size(M);
        if (theta~=0)
            M(l_out,:)=M(l_out,:)/M(l_out,1); % Copy row l_out, normalizing M(l_out,1) to 1
        end
        for k = 1:f % For all matrix rows
            if (k ~= l_out) % Other then l_out
                M(k,:)=M(k,:)-M(k,1)*M(l_out,:); % Set them equal to the original matrix Minus a multiple of normalized row l_out, making R(k,j_in)=0
            end
        end
        invB=M(1:mA,2:end);


        % Check if too many iterations are performed (increase n_max to 
        % allow more iterations)
        if(n == n_max)
            fprintf('Max number of iterations performed!\n\n');
            return
        end
    end 
end  

