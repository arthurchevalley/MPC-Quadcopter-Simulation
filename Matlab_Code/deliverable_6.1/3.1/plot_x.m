function = plot_x(sol,time)

figure
subplot(3,2,1)
hold on;grid on;
plot(time,sol.x(1,:),'Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Pitch speed",'Linewidth',1)

subplot(3,2,2)
hold on;grid on;
plot(time,sol.x(2,:),'Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("Pitch",'Linewidth',1)

subplot(3,2,3)
hold on;grid on;
plot(time,sol.x(3,:),'Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("X speed",'Linewidth',1)

subplot(3,2,4)
hold on;grid on;
plot(time,sol.x(1,:),'Linewidth',1)
xlabel("Time [s]",'Linewidth',1)
ylabel("X",'Linewidth',1)
end

