
function [F, L] = quad_discrete(time_step, quad)
%% Compute discrete-time dynamics 

import casadi.*

% Define variables 
alpha_dot = SX.sym('alpha_dot');
beta_dot = SX.sym('beta_dot');
gamma_dot = SX.sym('gamma_dot');
alpha = SX.sym('alpha');
beta = SX.sym('beta');
gamma = SX.sym('gamma');
x_dot = SX.sym('x_dot');
y_dot = SX.sym('y_dot');
z_dot = SX.sym('z_dot');
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');

state = [alpha_dot; beta_dot; gamma_dot; alpha; beta; gamma; x_dot; y_dot; z_dot; x; y; z];

u1 = SX.sym('u1');
u2 = SX.sym('u2');
u3 = SX.sym('u3');
u4 = SX.sym('u4');

u = [u1; u2; u3; u4];

ref_x = SX.sym('ref_x');
ref_y = SX.sym('ref_y');
ref_z = SX.sym('ref_z');
ref_gamma = SX.sym('ref_gamma');

ref = [ref_x; ref_y; ref_z; ref_gamma];

% Define system dynamics 

xdot = [ 0.5 * beta_dot * gamma_dot + 0.56*(u2-u4);...
        -0.5 * alpha_dot * gamma_dot + 0.56*(u3-u1);...
        11/15 * (u1 + u3 - u2 - u4);...
        alpha_dot;...
        beta_dot;...
        gamma_dot;...
        3.5*(u1 + u2 + u3 + u4) * sin(beta);...
        - 3.5*(u1 + u2 + u3 + u4) * sin(alpha) * cos(beta);...
        - 9.81 + 3.5*(u1 + u2 + u3 + u4) * cos(alpha) * cos(beta);...
        x_dot;...
        y_dot;...
        z_dot];

% Stage cost in order to minimize the distance to the reference and the
% speeds.
L = 40*(x - ref_x)^2         + ...     % minimize --> make x go to ref_x
    40*(y - ref_y)^2         + ...     % minimize --> make y go to ref_y
    200*(z - ref_z)^2        + ...     % minimize --> make z go to ref_z
    20*(gamma - ref_gamma)^2 + ...     % minimize --> make gamma go to ref_gamma
    1*u1^2                   + ...     % minimize first propeller speed
    1*u2^2                   + ...     % minimize second propeller speed
    1*u3^2                   + ...     % minimize third propeller speed
    1*u4^2                   + ...     % minimize fourth propeller speed
    10*alpha_dot^2           + ...     % minimize the angular speed 
    10*beta_dot^2            + ...     % minimize the angular speed
    10*gamma_dot^2           + ...     % minimize the angular speed
    10*x_dot^2               + ...     % minimize the linear speed
    10*y_dot^2               + ...     % minimize the linear speed
    50*z_dot^2 ;                       % minimize the linear speed
  

% Continous time dynamics 
f = Function('f', {state,u,ref}, {xdot, L});

% Formulate discrete time dynamics with Runge-Kutta 4 integrator
DT = time_step;
X0 = MX.sym('X0',12,1);
U = MX.sym('U',4,1);
REF = MX.sym('REF',4,1);
X = X0;
Q = 0;
[k1, k1_q] = f(X, U, REF);
[k2, k2_q] = f(X + DT/2 * k1, U, REF);
[k3, k3_q] = f(X + DT/2 * k2, U, REF);
[k4, k4_q] = f(X + DT * k3, U, REF);
X = X + DT/6 * (k1 + 2*k2 + 2*k3 + k4);
Q = Q + DT/6 * (k1_q + 2*k2_q + 2*k3_q + k4_q);

F = Function('F', {X0, U, REF}, {X});
L = Function('L', {X0, U, REF}, {Q});

end 