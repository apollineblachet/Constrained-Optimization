clear all 
close all
clc

%% Specification of Gradients + Specification of objective function

x0 = [-0.2;2.5]; % initial point
xl = [-5;-5]; % lower bounds
xu = [5;5]; % upper bounds

% -4x1 + 10 x2 >= 0 represented as A x <= b
A = [4 -10];
b = 0;
Aeq = zeros(0,2);
beq = zeros(0,1);

 % Parameters (in this case we do not use parameters)
p = [];

 % Call fmincon
options = optimoptions( 'fmincon',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,...
    'Display','none',...
    'Algorithm','interior-point',...
    'StepTolerance',1e-6); %tol on x

[sol,fval,exitflag,output]=fmincon( ...
    @objfungradHimmelblau, x0, ...
    A, b, Aeq, beq, ...
    xl, xu, ...
    @confungradHimmelblau1, ...
    options, ...
    p);

sol,fval,output


%% Plot Himmelblau problem with solution
x = -5:0.005:5;
y = -5:0.005:5;
[X,Y] = meshgrid(x,y);
F = (X.^2+Y-11).^2 + (X + Y.^2 - 7).^2;
v = [0:2:10 10:10:100 100:20:200];
[c,h]=contour(X,Y,F,v,'linewidth',2);
colorbar
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
    fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.2)
    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.2)
    plot(sol(1),sol(2),'r*')
hold off
xlim([-5,5])
ylim([-5,5])


%% Himmelblau function and non linear constraints function
function [f,dfdx] = objfungradHimmelblau(x,p)

tmp1 = x(1)*x(1)+x(2)-11;
tmp2 = x(1)+x(2)*x(2)-7;
f = tmp1*tmp1 + tmp2*tmp2;

% compute the gradient of f
if nargout > 1
    dfdx = zeros(2,1);
    dfdx(1,1) = 4*tmp1*x(1) + 2*tmp2;
    dfdx(2,1) = 2*tmp1 + 4*tmp2*x(2);
end
end

function [c,ceq,dcdx,dceqdx] = confungradHimmelblau1(x,p)

c = zeros(1,1);
ceq = zeros(0,1);

% Inequality constraints c(x) <= 0
tmp = x(1)+2;
c(1,1) = -(tmp*tmp - x(2));

% Compute constraint gradients
if nargout > 2
    dcdx = zeros(2,1);
    dceqdx = zeros(2,0);
    
    dcdx(1,1) = -2*tmp; % dc1dx1
    dcdx(2,1) = 1.0; % dc1dx2
end
end

