function [ctrl, traj] = My_ctrl_NMPC(quad) 

import casadi.*

% Input constraints

UMIN = 0;
UMAX = 1.5;

%% Formulate integrator (discrete-time model)

N = 100; % number of control intervals

% Time step for discretization
time_step = 0.26;

% Discretize the system
[F, L] = quad_discrete(time_step);

%% Define the optimization problem 

% N = ....; % MPC horizon [SET THIS VARIABLE]

opti = casadi.Opti(); % Optimization problem 

% Variables
X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4, N); % control trajectory 

% Parameters
X0 = opti.parameter(12,1); % initial state of the MPC problem (-> will change)
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]

% Dynamic constraints 
obj = 0;
for i = 1:N
    opti.subject_to(X(:,i+1) == F(X(:,i), U(:,i), REF));
    obj = obj + L(X(:,i), U(:,i), REF);
end 

% Input constraints 
opti.subject_to(UMIN <= U <= UMAX);

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

% Solve the optimization problem 
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g)); 
end




