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
%% Time vector init
for i = 1:Tf
    time = [time;(time(i)+Ts)];
end
%% x
sol.x(:,1) = [0;0;0;2];
for i = 1:Tf
    sol.ux(:,i) = mpc_x.get_u(sol.x(:,i));
    sol.x(:,i+1) = mpc_x.A*sol.x(:,i) + mpc_x.B*sol.ux(:,i);
end
sol.ux(:,i+1) = sol.ux(:,i);

%% y
sol.y(:,1) = [0;0;0;2];
for i = 1:Tf
    sol.uy(:,i) = mpc_y.get_u(sol.y(:,i));
    sol.y(:,i+1) = mpc_y.A*sol.y(:,i) + mpc_y.B*sol.uy(:,i);
end
sol.uy(:,i+1) = sol.uy(:,i);

%% z
sol.z(:,1) = [0;2];
for i = 1:Tf
    sol.uz(:,i) = mpc_z.get_u(sol.z(:,i));
    sol.z(:,i+1) = mpc_z.A*sol.z(:,i) + mpc_z.B*sol.uz(:,i);
end
sol.uz(:,i+1) = sol.uz(:,i);


%% yaw

sol.yaw(:,1) = [0;pi/4];
for i = 1:Tf
    sol.uyaw(:,i) = mpc_yaw.get_u(sol.yaw(:,i));
    sol.yaw(:,i+1) = mpc_yaw.A*sol.yaw(:,i) + mpc_yaw.B*sol.uyaw(:,i);
end

sol.uyaw(:,i+1) = sol.uyaw(:,i);

plot_evolution(sol,time);