% --------------------------
% Developer: Seung Jae Lee
% Date: 2/8/2021
% Abstract: This is the very first code to simulate dropping a ball
% Unit: cm, g, N
% --------------------------

clc
clear 
close all

% -------------------------------
% Global variables
% -------------------------------
grav_acc        = -981; % 981 cm/s^2

% user-input
sim_duration    = 0.2;  % unit: second 
dt              = 2e-4; % unit: second
sim_steps       = sim_duration/dt;
density         = 0.1;  % 0.1 g/cc 
% -------------------------------
% Define balls
% -------------------------------
b = struct();

% first ball
b(1).r  = 1;
b(1).m  = 4/3*pi*(b(1).r)^3*density;
b(1).cx = 5;    % center x
b(1).cz = 3;    % center z
b(1).vx = 0;    % velocity x
b(1).vz = 0;    % velocity z
b(1).fx = 0;    % force x
b(1).fz = 0;    % force z

% for theoretical solution
time = dt*[1:sim_steps];
theoretical_vz = b(1).vz + grav_acc*time;
theoretical_cz = b(1).cz + b(1).vz*time + 0.5*grav_acc*time.^2;

% -------------------------------
% Initial plot
% -------------------------------
% draw boundary
bndX = [0,0,10,10];
bndZ = [0,10,10,0];
patch(bndX,bndZ,'g','FaceAlpha',0.5);
axis equal 
grid on
hold on

% draw balls (circles)
Xcircle = cos([0:0.1:2*pi]);
Zcircle = sin([0:0.1:2*pi]);
GM(1) = patch(b(1).cx+Xcircle*b(1).r,b(1).cz+Zcircle*b(1).r,'b');
axis equal 

% -------------------------------
% Time integration (simulation)
% -------------------------------
array_vz = zeros(1,sim_steps);
array_cz = zeros(1,sim_steps);

for i = 1:sim_steps
    %--------------------------------------------
    % particle force update
    %--------------------------------------------    
    % global force update
    b(1).fz = b(1).m * grav_acc;

    %--------------------------------------------
    % particle motion update
    %--------------------------------------------
    % ball motion update
    b(1).vz     = b(1).vz + b(1).fz/b(1).m*dt;
    if (i==1) b(1).vz = b(1).vz/2; end % central time integration
    b(1).cz     = b(1).cz + b(1).vz*dt;
    array_vz(i) = b(1).vz;
    array_cz(i) = b(1).cz;
    
    %--------------------------------------------
    % update display 
    %--------------------------------------------
    set(GM(1),'XData',b(1).cx+Xcircle*b(1).r,'YData',b(1).cz+Zcircle*b(1).r)
    drawnow
end

% -------------------------------------
% Compare with the theoretical solution
% -------------------------------------
% velocity
figure
subplot(2,1,1)
hold on 
grid on 
plot(time,theoretical_vz,'r-')
plot(time,array_vz,'b--')
xlabel('Time (sec)')
ylabel('Velocity (cm/sec)')
legend('Theoretical','Simulation')

% displacement
subplot(2,1,2)
hold on 
grid on 
plot(time,theoretical_cz,'r-')
plot(time,array_cz,'b--')
xlabel('Time (sec)')
ylabel('Displacement (cm/sec)')
legend('Theoretical','Simulation')

sgtitle('Comparison with the theoretical solutions')




















