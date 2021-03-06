%--------------------------------------------
% Abstract: A super simple DEM code in 2D
% Unit: cm, g, N
%
%--------------------------------------------

clc
clear
close all

%--------------------------------------------
% Global variables
%--------------------------------------------
grav_acc        = -981;          % 981 cm/s^2

% user-input
sim_duration    = 0.5;         % unit: second
stiffness       = 1.3e5;       % 1.3e5 N/m = (kg*m/s^2)/m =(kg*cm/s^2)/cm 
density         = 0.1;         % 0.1 g/cm^3 
%--------------------------------------------
% Define balls (particles)
%--------------------------------------------
b = struct();

% first ball
b(1).r  = 1;
b(1).m  = 4/3*pi*b(1).r*b(1).r*b(1).r*density;
b(1).cx = 5;    % center x
b(1).cz = 5;    % center z
b(1).vx = 0;    % velocity x (change this for initial velocity)
b(1).vz = 0;    % velocity z (change this for initial velocity)
b(1).fx = 0;    % force x
b(1).fz = 0;    % force z

% second ball, third ball ... (likewise)


% compute the simulation step
min_r   = min(b(:).r);    % minimum ball size (radius)
dt      = 0.2 * sqrt(4/3*pi*min_r*min_r*min_r*density/stiffness); % time step size (unit: seconds)
sim_steps = ceil(sim_duration / dt);

%--------------------------------------------
% Initial plot
%--------------------------------------------
% draw boundary
bndX = [0,0,10,10]; 
bndZ = [0,10,10,0];
patch(bndX,bndZ,'g','FaceAlpha',0.5);
grid on;
axis square equal;
hold on;

% draw balls (circles)
Xcircle=cos([0:0.1:2*pi]);
Zcircle=sin([0:0.1:2*pi]);
GM(1) = patch (b(1).cx+b(1).r*Xcircle, b(1).cz+b(1).r*Zcircle,'b','FaceAlpha',1);    
%--------------------------------------------
% Time integration (simulation)
%--------------------------------------------
Nb = length(b);
time = dt*[1:sim_steps];
array_cz = zeros(1,length(1:sim_steps)); % array of b(1).cz

for i=1:sim_steps

    %--------------------------------------------
    % particle force update
    %--------------------------------------------
    % GLOBAL FORCE update
    b(1).fz = b(1).m * grav_acc; 
    
    % CONTACT FORCE update  
    % Check boundaries 
    % Assumption: the boundary is rectangle. If not rectangle, a more
    %             sophisticated method needs to be implemented.
    urn = (b(1).cz-b(1).r) - min(bndZ);
    if urn < 0
        cnt_force = stiffness * urn;
        b(1).fz = b(1).fz - cnt_force;
    end
    
    %--------------------------------------------
    % particle motion update
    %--------------------------------------------
    b(1).vz = b(1).vz + b(1).fz/b(1).m*dt;    
    b(1).cz = b(1).cz + b(1).vz*dt;
    array_cz(i) = b(1).cz;
    
    %--------------------------------------------
    % update display 
    %--------------------------------------------
    set(GM(1),'XData',b(1).cx+b(1).r*Xcircle, 'YData', b(1).cz+b(1).r*Zcircle);
    drawnow
    
end

%--------------------------------------------
% Plot displacement history
%--------------------------------------------
figure;  
hold on
grid on
plot(time,array_cz,'b')
title('Simulated Displacement')
xlabel('Time (sec)')
ylabel('Displacement (cm)')

