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
sim_duration    = 0.3;  % unit: second 
stiffness       = 1.3e5;% unit: N/m
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
b(1).vx = 150;    % velocity x
b(1).vz = 200;    % velocity z
b(1).fx = 0;    % force x
b(1).fz = 0;    % force z

% compute the simulation step
min_r   =       min(b(:).r);    % minimum ball size (radius)
mass_min=       4/3*pi*min_r*min_r*min_r*density;
dt      =       0.2*sqrt(mass_min/stiffness); % unit: second
sim_steps       = round(sim_duration/dt);
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
axis ([-2 12 -2 12])
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

array_cz = zeros(1,sim_steps);

for i = 1:sim_steps
    %--------------------------------------------
    % particle force update
    %--------------------------------------------    
    % global force update
    b(1).fz = b(1).m * grav_acc;
    b(1).fx = 0;

    % CONTACT FORCE update  
    % Check boundaries 
    % Assumption: the boundary is rectangle. If not rectangle, a more
    %             sophisticated method needs to be implemented.    
    % bottom boundary
    if (b(1).cz - b(1).r) < min(bndZ)
        cnt_force = stiffness * (min(bndZ) - (b(1).cz - b(1).r));
        b(1).fz = b(1).fz + cnt_force;
    end
    
    % top boundary
    if (b(1).cz + b(1).r) > max(bndZ)
        cnt_force = stiffness * (max(bndZ) - (b(1).cz + b(1).r));
        b(1).fz = b(1).fz + cnt_force;
    end    
    
    % left boundary
    if (b(1).cx - b(1).r) < min(bndX)
        cnt_force = stiffness * (min(bndX) - (b(1).cx - b(1).r));
        b(1).fx = b(1).fx + cnt_force;
    end

    % right boundary
    if (b(1).cx + b(1).r) > max(bndX)
        cnt_force = stiffness * (max(bndX) - (b(1).cx + b(1).r));
        b(1).fx = b(1).fx + cnt_force;
    end    
    
    %--------------------------------------------
    % particle motion update
    %--------------------------------------------
    % ball motion update
    b(1).vx     = b(1).vx + b(1).fx/b(1).m*dt;    
    b(1).vz     = b(1).vz + b(1).fz/b(1).m*dt;    
    b(1).cx     = b(1).cx + b(1).vx*dt;
    b(1).cz     = b(1).cz + b(1).vz*dt;
    array_cz(i) = b(1).cz;
    
    %--------------------------------------------
    % update display 
    %--------------------------------------------
    set(GM(1),'XData',b(1).cx+Xcircle*b(1).r,'YData',b(1).cz+Zcircle*b(1).r)
    drawnow
end

%--------------------------------------------
% Plot displacement history
%--------------------------------------------

% displacement
figure
grid on 
plot(time,array_cz,'b-')
xlabel('Time (sec)')
ylabel('Displacement (cm)')






















