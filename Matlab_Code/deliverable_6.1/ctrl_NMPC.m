function [ctrl, traj] = ctrl_NMPC(quad) 

import casadi.*

% Input constraints
UMIN = zeros(4,1);
UMAX = 1.5*ones(4,1);

N = 30; % MPC horizon 

% Time step for discretization
h = 0.26;

% Dicretization using Runga-Kutta 4
f_discrete = @(x,u) RK4(x,u,h,quad);

%% Define the optimization problem 

opti = casadi.Opti(); % Optimization problem 

% Variables
X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4, N); % control trajectory 

% Parameters
X0 = opti.parameter(12,1); % initial state of the MPC problem (-> will change)
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]

obj = 0;
for k=1:N % loop over control intervals
  
  % Dynamics constraints 
  opti.subject_to(X(:,k+1) == f_discrete(X(:,k), U(:,k)));
  
  % Input constraints 
  opti.subject_to(UMIN <= U(:,k) <= UMAX);
  
  % Objective of the optimization
  obj = obj + U(:,k)'*U(:,k);                           % Minimize the propeller speeds (u1, u2, u3, u4)
  obj = obj + 5*(X(10,k) - REF(1))'*(X(10,k) - REF(1)); % Minimize --> make state x go to ref_x
  obj = obj + 5*(X(11,k) - REF(2))'*(X(11,k) - REF(2)); % Minimize --> make state y go to ref_y 
  obj = obj + 5*(X(12,k) - REF(3))'*(X(12,k) - REF(3)); % Minimize --> make state z go to ref_z
  obj = obj + 5*( X(6,k) - REF(4))'*( X(6,k) - REF(4)); % Minimize --> make state yaw go to ref_yaw
end

% Initial conditions 
opti.subject_to(X(:,1) == X0);

% Objective
opti.minimize(obj);

ctrl = @(x,ref) eval_ctrl(x, ref, opti, X0, REF, X, U);
end

function u = eval_ctrl(x, ref, opti, X0, REF, X, U) 

% Set the initial state and reference 
opti.set_value(X0, x);
opti.set_value(REF, ref);

% ipopt: open source interior point solver for nn-linear problems
ops = struct('ipopt', struct('print_level',0, 'tol', 1e-3), 'print_time', false); 
% ops = struct('ipopt', struct('tol', 1e-3), 'print time', false); 
opti.solver('ipopt', ops);

% Solve the optimization problem giving optimal input
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g)); 
end








