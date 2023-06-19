function [ctrl, traj] = ctrl NMPC(quad)
import casadi.*
opti = casadi.Opti(); % Optimization problem
N = 20 % MPC horizon [SET THIS VARIABLE]
% −−−− decision variables −−−−−−−−−
X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4, N); % control trajectory (throttle, brake)
X0 = opti.parameter(12,1); % initial state
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% YOUR CODE HERE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
ctrl = @(x,ref) eval ctrl(x, ref, opti, X0, REF, X, U);
end
function u = eval ctrl(x, ref, opti, X0, REF, X, U)
% −−−− Set the initial state and reference −−−−
opti.set value(X0, x);
opti.set value(REF, ref);
% −−−− Setup solver NLP −−−−−−
ops = struct('ipopt', struct('print level',0, 'tol', 1e−3), 'print time', false);
opti.solver('ipopt', ops);
% −−−− Solve the optimization problem −−−−
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');
u = opti.value(U(:,1));
% Use the current solution to speed up the next optimization
opti.set initial(sol.value variables());
opti.set initial(opti.lam g, sol.value(opti.lam g));
end