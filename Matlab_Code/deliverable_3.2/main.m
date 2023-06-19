
%--------------------------------------------------------------------------
%% ----------------------- MAIN of deliverable 3.2 ------------------------
% -------------------------------------------------------------------------

clc
clear
Ts = 1/5;
quad = Quad(Ts);
[xs,us] = quad.trim();
sys = quad.linearize(xs, us); % Linearize the nonlinear model
sys_transformed = sys * inv(quad.T); % New system is A * x + B * inv(T) * v
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);


%% ------------------ sys_x ---------------------------------
mpc_x = MPC_Control_x(sys_x, Ts);

x_position_reference = -2;
x0_sys_x = [0 0 0 0]';
sol_x.x(:,1) = x0_sys_x;
i = 1;

while norm(sol_x.x(4,end) - x_position_reference) > 1e-3 
    
    sol_x.u(:,i) = mpc_x.get_u(sol_x.x(:,i), x_position_reference); 
    sol_x.x(:,i+1) = mpc_x.A * sol_x.x(:,i) + mpc_x.B * sol_x.u(:,i);
    i = i + 1
end 

time_x = 1:length(sol_x.x);
figure()
subplot(2,2,1)
plot(time_x*Ts, sol_x.x(1,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('Pitch angular velocity in [rad/s]')

subplot(2,2,2)
plot(time_x*Ts, sol_x.x(2,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('Pitch angle in [rad]')

subplot(2,2,3)
plot(time_x*Ts, sol_x.x(3,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('linear velocity along x in [m/s]')

subplot(2,2,4)
plot(time_x*Ts, sol_x.x(4,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('position along x in [m]')

sgt = sgtitle('Evolution of the 4 states of sys_x','Color','black');
sgt.FontSize = 20;

figure()
time_u = 1:length(sol_x.u);
plot(time_u*Ts, sol_x.u,'LineWidth',2)    
xlabel('time in [s]')
ylabel('Pitch body moment')
sgt = sgtitle('Evolution of the input of sys_x','Color','black');
sgt.FontSize = 20;  

%% ---------------------------- sys_y ---------------------------------

mpc_y = MPC_Control_y(sys_y, Ts);

y_position_reference = -2;
x0_sys_y = [0 0 0 0]';

sol_y.x(:,1) = x0_sys_y;
i = 1;

while norm(sol_y.x(4,end) - y_position_reference) > 1e-3 
    sol_y.u(:,i) = mpc_y.get_u(sol_y.x(:,i), y_position_reference);
    sol_y.x(:,i+1) = mpc_y.A * sol_y.x(:,i) + mpc_y.B * sol_y.u(:,i);
    i = i + 1
end 

figure()
time_x = 1:length(sol_y.x);
subplot(2,2,1)
plot(time_x*Ts, sol_y.x(1,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('Roll angular velocity in [rad/s]')

subplot(2,2,2)
plot(time_x*Ts, sol_y.x(2,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('Roll angle in [rad]')

subplot(2,2,3)
plot(time_x*Ts, sol_y.x(3,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('linear velocity along y in [m/s]')

subplot(2,2,4)
plot(time_x*Ts, sol_y.x(4,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('position along y in [m]')

sgt = sgtitle('Evolution of the 4 states of sys_y','Color','black');
sgt.FontSize = 20;

figure()
time_u = 1:length(sol_y.u);
plot(time_u*Ts, sol_y.u,'LineWidth',2)    
xlabel('time in [s]')
ylabel('Roll body moment')
sgt = sgtitle('Evolution of the input of sys_y','Color','black');
sgt.FontSize = 20;  


%% ---------------------------- sys_z ---------------------------------

mpc_z = MPC_Control_z(sys_z, Ts);
x0_sys_z = [0 0]';

z_position_reference = -2;
sol_z.x(:,1) = x0_sys_z;
i = 1;

% pour 1e-2, j'ai 40/5 = 8 sec 
% pour 1e-3, j'ai 51/5 = 10 sec 
while norm(sol_z.x(2,end) - z_position_reference) > 1e-3 
    sol_z.u(:,i) = mpc_z.get_u(sol_z.x(:,i), z_position_reference);
    sol_z.x(:,i+1) = mpc_z.A * sol_z.x(:,i) + mpc_z.B * sol_z.u(:,i);
    i = i + 1
end 

figure()
time_x = 1:length(sol_z.x);
subplot(1,2,1)
plot(time_x*Ts, sol_z.x(1,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('linear velocity along z in [m/s]')

subplot(1,2,2)
plot(time_x*Ts, sol_z.x(2,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('position along z in [m]')

sgt = sgtitle('Evolution of the 2 states of sys_z','Color','black');
sgt.FontSize = 20;

figure()
time_u = 1:length(sol_z.u);
plot(time_u*Ts, sol_z.u,'LineWidth',2)    
xlabel('time in [s]')
ylabel('Force in [N]')
sgt = sgtitle('Evolution of the input of sys_z','Color','black');
sgt.FontSize = 20;

%% ---------------------------- sys_yaw ---------------------------------

mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);
x0_sys_yaw = [0 0]';

yaw_position_reference = pi/4;

sol_yaw.x(:,1) = x0_sys_yaw;
i = 1;

while norm(sol_yaw.x(2,end) - yaw_position_reference) > 1e-3 
    sol_yaw.u(:,i) = mpc_yaw.get_u(sol_yaw.x(:,i), yaw_position_reference);
    sol_yaw.x(:,i+1) = mpc_yaw.A * sol_yaw.x(:,i) + mpc_yaw.B * sol_yaw.u(:,i);
    i = i + 1
end 

figure()
time_x = 1:length(sol_yaw.x);
subplot(1,2,1)
plot(time_x*Ts, sol_yaw.x(1,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('Yaw angular velocity in [rad/s]')

subplot(1,2,2)
plot(time_x*Ts, sol_yaw.x(2,:),'LineWidth',2) 
xlabel('time in [s]')
ylabel('Yaw angle in [rad]')

sgt = sgtitle('Evolution of the 2 states of sys_yaw','Color','black');
sgt.FontSize = 20;

figure()
time_u = 1:length(sol_yaw.u);
plot(time_u*Ts, sol_yaw.u,'LineWidth',2)    
xlabel('time in [s]')
ylabel('Yaw body moment')
sgt = sgtitle('Evolution of the input of sys_yaw','Color','black');
sgt.FontSize = 20; 

