function plot_evolution(sol,time)

% X plots
figure('Name','X subsystem','NumberTitle','off')
subplot(3,2,1)
grid on;
plot(time,sol.x(1,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Pitch speed [rad/s]",'Linewidth',1)
subplot(3,2,2)
grid on;
plot(time,sol.x(2,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Pitch [rad",'Linewidth',1)
subplot(3,2,3)
grid on;
plot(time,sol.x(3,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("X speed [m/s]",'Linewidth',1)
subplot(3,2,4)
grid on;
plot(time,sol.x(4,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("X [m]",'Linewidth',1)
grid on;
subplot(3,2,5)
plot(time,sol.ux,'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("M_\alpha [Nm]",'Linewidth',1)

% y plots
figure('Name','Y subsystem','NumberTitle','off')
subplot(3,2,1)
grid on;
plot(time,sol.y(1,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Roll velocity [rad/s]",'Linewidth',1)
subplot(3,2,2)
grid on;
plot(time,sol.y(2,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Roll [rad]",'Linewidth',1)
subplot(3,2,3)
grid on;
plot(time,sol.y(3,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Y speed [m/s]",'Linewidth',1)
subplot(3,2,4)
grid on;
plot(time,sol.y(4,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Y [m]",'Linewidth',1)
subplot(3,2,5)
grid on;
plot(time,sol.uy,'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("M_\beta [Nm]",'Linewidth',1)

%z
figure('Name','Z subsystem','NumberTitle','off')
subplot(3,2,1)
grid on;
plot(time,sol.z(1,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Z speed [m/s]",'Linewidth',1)
subplot(3,2,2)
grid on;
plot(time,sol.z(2,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Z [m]",'Linewidth',1)
subplot(3,2,3)
grid on;
plot(time,sol.uz,'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("F [N]",'Linewidth',1)


% yaw
figure('Name','Yaw subsystem','NumberTitle','off')
subplot(3,2,1)
grid on;
plot(time,sol.yaw(1,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Yaw velocity [rad/s]",'Linewidth',1)
subplot(3,2,2)
grid on;
plot(time,sol.yaw(2,:),'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Yaw [rad]",'Linewidth',1)
subplot(3,2,3)
grid on;
plot(time,sol.uyaw,'r','Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("M_\gamma [Nm]",'Linewidth',1)
end

