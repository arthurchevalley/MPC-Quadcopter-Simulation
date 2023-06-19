function [ctrl, traj] = ctrl_MPC(quad)

opti = casadi.Opti(); % Optimization problem

N = 30; % MPC horizon 
% −−−− decision variables −−−−−−−−−

X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4, N); % control trajectory (throttle, brake)
X0 = opti.parameter(12,1); % initial state
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]

% Time step in order to discretize the system using RK4
h = 0.2;

% Weights matrices for the objective function
Q = 5*eye(12);
R = eye(4);

% Discretization
discrete = @(x,u) RK4(x,u,h,quad);

% Reference vector in the same size as the state vector
r = [zeros(5,1); REF(4,1); zeros(3,1); REF(1:3,1)];

% Dynamic constraints 
cost = (X(:,1)-r)'*Q*(X(:,1)-r)+U(:,1)'*R*U(:,1);
opti.subject_to(U(:,1) <= 1.5*ones(4,1));
opti.subject_to(U(:,1) >= zeros(4,1));
opti.subject_to(X(:,2) == discrete(X(:,1),U(:,1)));
opti.subject_to(X(:,1)== X0);

for j=2:N
    opti.subject_to(U(:,j) <= 1.5*ones(4,1));
    opti.subject_to(U(:,j) >= zeros(4,1));
    opti.subject_to(X(:,j+1) == discrete(X(:,j),U(:,j)));
    cost = cost + (X(:,j)-r)'*Q*(X(:,j)-r)+U(:,j)'*R*U(:,j);
end
opti.minimize(cost);

ctrl = @(x,ref) eval_ctrl(x, ref, opti, X0, REF, X, U);
end
function u = eval_ctrl(x, ref, opti, X0, REF, X, U)
% −−−− Set the initial state and reference −−−−
opti.set_value(X0, x);
opti.set_value(REF, ref);
% −−−− Setup solver NLP −−−−−−
ops = struct('ipopt', struct('print_level',0, 'tol', 1/1000), 'print_time', false);
opti.solver('ipopt', ops);
% −−−− Solve the optimization problem −−−−
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');
u = opti.value(U(:,1));
% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end