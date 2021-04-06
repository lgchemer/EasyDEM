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
sim_duration    = 1;  % unit: second 
stiffness       = 1.3e5;% unit: N/m
density         = 0.1;  % 0.1 g/cc
cnt_dmp_ratio   = 0.1;    % This cannot exceed 1 (= Critical damping)
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

% second ball
b(2).r  = 1;
b(2).m  = 4/3*pi*(b(2).r)^3*density;
b(2).cx = 2;    % center x
b(2).cz = 3;    % center z
b(2).vx = 100;    % velocity x
b(2).vz = 0;    % velocity z
b(2).fx = 0;    % force x
b(2).fz = 0;    % force z

% compute the simulation step
min_r   =       min(b(:).r);                    % minimum ball size (radius)
mass_min=       4/3*pi*min_r*min_r*min_r*density;
dt      =       0.2*sqrt(mass_min/stiffness);   % unit: second
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
GM(2) = patch(b(2).cx+Xcircle*b(2).r,b(2).cz+Zcircle*b(2).r,'r');
axis equal 

% -------------------------------
% Time integration (simulation)
% -------------------------------

array_cz = zeros(2,sim_steps);

for i = 1:sim_steps
    %--------------------------------------------
    % particle force update
    %--------------------------------------------    
    % global force update
    b(1).fz = b(1).m * grav_acc;
    b(1).fx = 0;
    
    b(2).fz = b(2).m * grav_acc;
    b(2).fx = 0;
    
    % CONTACT FORCE update  
    % Check boundaries 
    % Assumption: the boundary is rectangle. If not rectangle, a more
    %             sophisticated method needs to be implemented.    
    %--------------------------------------------
    % Ball 1 interaction with boundary
    %--------------------------------------------      
    % bottom boundary
    if (b(1).cz - b(1).r) < min(bndZ)
        cnt_force = stiffness * (min(bndZ) - (b(1).cz - b(1).r));
        b(1).fz = b(1).fz + cnt_force;  % add contact force (Fne)
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vz;
        b(1).fz =  b(1).fz + cnt_damp;  % add contact damping (Fnd)
    end
    
    % top boundary
    if (b(1).cz + b(1).r) > max(bndZ)
        cnt_force = stiffness * (max(bndZ) - (b(1).cz + b(1).r));
        b(1).fz = b(1).fz + cnt_force;  % add contact force (Fne)
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vz; 
        b(1).fz =  b(1).fz + cnt_damp;  % add contact damping (Fnd)
    end    
    
    % left boundary
    if (b(1).cx - b(1).r) < min(bndX)
        cnt_force = stiffness * (min(bndX) - (b(1).cx - b(1).r));
        b(1).fx = b(1).fx + cnt_force;  % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vx; 
        b(1).fx = b(1).fx + cnt_damp;  % add contact damping (Fnd)        
    end

    % right boundary
    if (b(1).cx + b(1).r) > max(bndX)
        cnt_force = stiffness * (max(bndX) - (b(1).cx + b(1).r));
        b(1).fx = b(1).fx + cnt_force;
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vx; 
        b(1).fx = b(1).fx + cnt_damp;  % add contact damping (Fnd)           
    end    
    
    %--------------------------------------------
    % Ball 2 interaction with boundary
    %--------------------------------------------       
    if (b(2).cz - b(2).r) < min(bndZ)
        cnt_force = stiffness * (min(bndZ) - (b(2).cz - b(2).r));
        b(2).fz = b(2).fz + cnt_force;  % add contact force (Fne)
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vz;
        b(2).fz =  b(2).fz + cnt_damp;  % add contact damping (Fnd)
    end
    
    % top boundary
    if (b(2).cz + b(2).r) > max(bndZ)
        cnt_force = stiffness * (max(bndZ) - (b(2).cz + b(2).r));
        b(2).fz = b(2).fz + cnt_force;  % add contact force (Fne)
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vz; 
        b(2).fz =  b(2).fz + cnt_damp;  % add contact damping (Fnd)
    end    
    
    % left boundary
    if (b(2).cx - b(2).r) < min(bndX)
        cnt_force = stiffness * (min(bndX) - (b(2).cx - b(2).r));
        b(2).fx = b(2).fx + cnt_force;  % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vx; 
        b(2).fx = b(2).fx + cnt_damp;  % add contact damping (Fnd)        
    end

    % right boundary
    if (b(2).cx + b(2).r) > max(bndX)
        cnt_force = stiffness * (max(bndX) - (b(2).cx + b(2).r));
        b(2).fx = b(2).fx + cnt_force;
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vx; 
        b(2).fx = b(2).fx + cnt_damp;  % add contact damping (Fnd)           
    end 
    
    %--------------------------------------------
    % Two balls interacts each other
    %--------------------------------------------      
    cnt_n = [b(2).cx-b(1).cx b(2).cz-b(1).cz];      % contact normal from ball 1 to 2
    urn   = vecnorm(cnt_n) - b(2).r - b(1).r;       % if some overlap, urn < 0
    
    if urn < 0
        cnt_force   = stiffness * urn;                      % contact force
        b(1).fx     = b(1).fx + cnt_force * cnt_n(1);
        b(1).fz     = b(1).fz + cnt_force * cnt_n(2);
        b(2).fx     = b(2).fx + cnt_force * cnt_n(1) * -1;
        b(2).fz     = b(2).fz + cnt_force * cnt_n(2) * -1;  
        
        crit_dmp    = 2 * sqrt(b(1).m*b(2).m*stiffness/(b(1).m+b(2).m));
        cnt_damp    = cnt_dmp_ratio * crit_dmp * dot([b(2).vx-b(1).vx b(2).vz-b(1).vz],cnt_n);
        b(1).fx     = b(1).fx + cnt_damp * cnt_n(1);
        b(1).fz     = b(1).fz + cnt_damp * cnt_n(2);
        b(2).fx     = b(2).fx + cnt_damp * cnt_n(1) * -1;
        b(2).fz     = b(2).fz + cnt_damp * cnt_n(2) * -1;          
    end
    
    %--------------------------------------------
    % particle motion update
    %--------------------------------------------
    % ball motion update
    b(1).vx     = b(1).vx + b(1).fx/b(1).m*dt;    
    b(1).vz     = b(1).vz + b(1).fz/b(1).m*dt;    
    b(1).cx     = b(1).cx + b(1).vx*dt;
    b(1).cz     = b(1).cz + b(1).vz*dt;
    array_cz(1,i) = b(1).cz;
    
    b(2).vx     = b(2).vx + b(2).fx/b(2).m*dt;    
    b(2).vz     = b(2).vz + b(2).fz/b(2).m*dt;    
    b(2).cx     = b(2).cx + b(2).vx*dt;
    b(2).cz     = b(2).cz + b(2).vz*dt;
    array_cz(2,i) = b(2).cz;
    
    
    %--------------------------------------------
    % update display 
    %--------------------------------------------
    set(GM(1),'XData',b(1).cx+Xcircle*b(1).r,'YData',b(1).cz+Zcircle*b(1).r)
    set(GM(2),'XData',b(2).cx+Xcircle*b(2).r,'YData',b(2).cz+Zcircle*b(2).r)
    title(gca, ['Time = ', num2str(dt*i),' seconds']);
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






















