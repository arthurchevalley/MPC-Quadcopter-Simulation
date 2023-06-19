clc
Tf = 100; % simulation time
Ts = 1/5;
time = [0];
quad = Quad(Ts);
%Tf = 1.0;                          %time of simulation
%x0 = zeros(12,1);                  %start state
%u = [0;0;1.5;1.5];                 %inputs
%sim = ode45(@(t,x) quad.f(x,u),[0,Tf],x0);
%quad.plot(sim,u)

[xs,us] = quad.trim();
sys = quad.linearize(xs,us);

[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);
mpc_x = MPC_Control_x(sys_x,Ts);
mpc_y = MPC_Control_y(sys_y, Ts);
mpc_z = MPC_Control_z(sys_z, Ts);
mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);

%% ref 
sim = quad.sim(mpc_x,mpc_y,mpc_z, mpc_yaw);
quad.plot(sim)

