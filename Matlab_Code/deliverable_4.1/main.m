
%--------------------------------------------------------------------------
%% ----------------------- MAIN of deliverable 4.1 ------------------------
% -------------------------------------------------------------------------
% Team members: - Pereira Portela Tifanny 
%               - Chevalley Arthur 
%               - Mermoud Paco
%
% Date: Autumn 2020
%
%
% The aim of this deliverable is to test the controllers in a simulation
% tracking a given path
%
%
% -------------------------------------------------------------------------

clc
clear
Ts = 1/5;
quad = Quad(Ts);
[xs, us] = quad.trim();
sys = quad.linearize(xs, us); % Linearize the nonlinear model
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);

% Subsystem controller initialization
mpc_x = MPC_Control_x(sys_x, Ts);
mpc_y = MPC_Control_y(sys_y, Ts);
mpc_z = MPC_Control_z(sys_z, Ts);
mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);

% Simulate the system in order to follow the path
sim = quad.sim(mpc_x, mpc_y, mpc_z, mpc_yaw); 
quad.plot(sim);