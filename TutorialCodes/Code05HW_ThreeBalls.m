%--------------------------------------------
% Developer: Seung Jae Lee
% Institution: Florida Interantional University
% Date: 3/1/2021
% Abstract: A super simple DEM code in 2D
%           This is the very first code to simulate dropping a ball
% Unit: cm, g, N
%--------------------------------------------

clc
clear
close all

%--------------------------------------------
% Global variables
%--------------------------------------------
grav_acc        = -981;          % 981 cm/s^2

% user-input
sim_duration    = 0.3;           % unit: second
stiffness       = 1.3e5;         % 1.3e5 N/m = (kg*m/s^2)/m =(kg*cm/s^2)/cm 
density         = 0.1;           % 0.1 g/cm^3 
cnt_dmp_ratio   = 0.1;           % This cannot exceed 1 (= Critical damping)
%--------------------------------------------
% Define balls (particles)
%--------------------------------------------
b = struct();

% first ball
b(1).r  = 1;
b(1).m  = 4/3*pi*b(1).r*b(1).r*b(1).r*density;
b(1).cx = 5;    % center x
b(1).cz = 3;    % center z
b(1).vx = 0;    % velocity x (change this for initial velocity)
b(1).vz = 50;    % velocity z (change this for initial velocity)
b(1).fx = 0;    % force x
b(1).fz = 0;    % force z

% second ball, third ball ... (likewise)
b(2).r  = 1;
b(2).m  = 4/3*pi*b(1).r*b(1).r*b(1).r*density;
b(2).cx = 2;  % center x
b(2).cz = 3;    % center z
b(2).vx = 100;    % velocity x (change this for initial velocity)
b(2).vz = 0;    % velocity z (change this for initial velocity)
b(2).fx = 0;    % force x
b(2).fz = 0;    % force z

% third ball, third ball ... (likewise)
b(3).r  = 1;
b(3).m  = 4/3*pi*b(1).r*b(1).r*b(1).r*density;
b(3).cx = 2;  % center x
b(3).cz = 6;    % center z
b(3).vx = 0;    % velocity x (change this for initial velocity)
b(3).vz = 0;    % velocity z (change this for initial velocity)
b(3).fx = 0;    % force x
b(3).fz = 0;    % force z

% compute the simulation step
min_r   = min(struct2table(b).r);    % minimum ball size (radius)
dt      = 0.1 * sqrt(4/3*pi*min_r*min_r*min_r*density/stiffness); % time step size (unit: seconds)
sim_steps = ceil(sim_duration / dt);

%--------------------------------------------
% Initial plot
%--------------------------------------------
% draw boundary
bndX    = [0,0,10,10]; 
bndZ    = [0,10,10,0];
patch(bndX,bndZ,'g','FaceAlpha',0.5);
grid on;
axis ([-2 12 -2 12])
axis square equal;
hold on;

% draw balls (circles)
Xcircle=cos([0:0.1:2*pi]);
Zcircle=sin([0:0.1:2*pi]);
GM(1) = patch (b(1).cx+b(1).r*Xcircle, b(1).cz+b(1).r*Zcircle,'b','FaceAlpha',1);    
GM(2) = patch (b(2).cx+b(2).r*Xcircle, b(2).cz+b(2).r*Zcircle,'r','FaceAlpha',1);
GM(3) = patch (b(3).cx+b(3).r*Xcircle, b(3).cz+b(3).r*Zcircle,'g','FaceAlpha',1);

%--------------------------------------------
% Time integration (simulation)
%--------------------------------------------
Nb = length(b);
time = dt*[1:sim_steps];
array_cz = zeros(3,length(1:sim_steps)); % array of b(#).cz

for i=1:sim_steps

    %--------------------------------------------
    % particle force update
    %--------------------------------------------
    % GLOBAL FORCE update
    b(1).fz = b(1).m * grav_acc; 
    b(1).fx = 0;

    b(2).fz = b(2).m * grav_acc; 
    b(2).fx = 0;    
    
    b(3).fz = b(3).m * grav_acc; 
    b(3).fx = 0;    
        
    % CONTACT FORCE update  
    % Check boundaries 
    % Assumption: the boundary is rectangle. If not rectangle, a more
    %             sophisticated method needs to be implemented.

    %--------------------------------------------
    % Ball 1 interaction with boundary
    %--------------------------------------------    
    % bottom boundary
    if (b(1).cz-b(1).r) < min(bndZ)     
        cnt_force = stiffness * (min(bndZ) - (b(1).cz-b(1).r)); 
        b(1).fz = b(1).fz + cnt_force;  % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vz ; 
        b(1).fz = b(1).fz + cnt_damp;  % add damping force
    end
    
    % top boundary
    if (b(1).cz+b(1).r) > max(bndZ)     
        cnt_force = stiffness * (max(bndZ) - (b(1).cz+b(1).r));
        b(1).fz = b(1).fz + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vz ; 
        b(1).fz = b(1).fz + cnt_damp;  % add damping force        
    end

    % left boundary
    if (b(1).cx-b(1).r) < min(bndX)     
        cnt_force = stiffness * (min(bndX) - (b(1).cx-b(1).r));
        b(1).fx = b(1).fx + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vx ; 
        b(1).fx = b(1).fx + cnt_damp;  % add damping force               
    end     
    
    % right boundary
    if (b(1).cx+b(1).r) > max(bndX)     
        cnt_force = stiffness * (max(bndX) - (b(1).cx+b(1).r));
        b(1).fx = b(1).fx + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(1).m) * -1 * b(1).vx ; 
        b(1).fx = b(1).fx + cnt_damp;  % add damping force               
    end    

    %--------------------------------------------
    % Ball 2 interaction with boundary
    %--------------------------------------------     
    % bottom boundary
    if (b(2).cz-b(2).r) < min(bndZ)     
        cnt_force = stiffness * (min(bndZ) - (b(2).cz-b(2).r)); 
        b(2).fz = b(2).fz + cnt_force;  % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vz ; 
        b(2).fz = b(2).fz + cnt_damp;  % add contact damping 
    end
    
    % top boundary
    if (b(2).cz+b(2).r) > max(bndZ)     
        cnt_force = stiffness * (max(bndZ) - (b(2).cz+b(2).r));
        b(2).fz = b(2).fz + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vz ; 
        b(2).fz = b(2).fz + cnt_damp; % add contact damping         
    end

    % left boundary
    if (b(2).cx-b(2).r) < min(bndX)     
        cnt_force = stiffness * (min(bndX) - (b(2).cx-b(2).r));
        b(2).fx = b(2).fx + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vx ; 
        b(2).fx = b(2).fx + cnt_damp; % add contact damping                
    end     
    
    % right boundary
    if (b(2).cx+b(2).r) > max(bndX)     
        cnt_force = stiffness * (max(bndX) - (b(2).cx+b(2).r));
        b(2).fx = b(2).fx + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(2).m) * -1 * b(2).vx ; 
        b(2).fx = b(2).fx + cnt_damp; % add contact damping                
    end        

    %--------------------------------------------
    % Ball 3 interaction with boundary
    %--------------------------------------------     
    % bottom boundary
    if (b(3).cz-b(3).r) < min(bndZ)     
        cnt_force = stiffness * (min(bndZ) - (b(3).cz-b(3).r)); 
        b(3).fz = b(3).fz + cnt_force;  % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(3).m) * -1 * b(3).vz ; 
        b(3).fz = b(3).fz + cnt_damp;  % add contact damping 
    end
    
    % top boundary
    if (b(3).cz+b(3).r) > max(bndZ)     
        cnt_force = stiffness * (max(bndZ) - (b(3).cz+b(3).r));
        b(3).fz = b(3).fz + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(3).m) * -1 * b(3).vz ; 
        b(3).fz = b(3).fz + cnt_damp; % add contact damping         
    end

    % left boundary
    if (b(3).cx-b(3).r) < min(bndX)     
        cnt_force = stiffness * (min(bndX) - (b(3).cx-b(3).r));
        b(3).fx = b(3).fx + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(3).m) * -1 * b(3).vx ; 
        b(3).fx = b(3).fx + cnt_damp; % add contact damping                
    end     
    
    % right boundary
    if (b(3).cx+b(3).r) > max(bndX)     
        cnt_force = stiffness * (max(bndX) - (b(3).cx+b(3).r));
        b(3).fx = b(3).fx + cnt_force; % add contact force
        cnt_damp = cnt_dmp_ratio * 2 * sqrt(stiffness*b(3).m) * -1 * b(3).vx ; 
        b(3).fx = b(3).fx + cnt_damp; % add contact damping                
    end         
    
    %--------------------------------------------
    % Two balls interacts each other
    %--------------------------------------------       
    cnt_n   = [b(2).cx-b(1).cx b(2).cz-b(1).cz];    % contact normal from ball 1 to 2
    urn     = vecnorm(cnt_n) - b(2).r - b(1).r;     % if some overlap, urn < 0
    if urn < 0  % Balls in contact
        cnt_force   = stiffness * urn;                      % contact force
        b(1).fx     = b(1).fx + cnt_force * cnt_n(1);
        b(1).fz     = b(1).fz + cnt_force * cnt_n(2);
        b(2).fx     = b(2).fx + cnt_force * cnt_n(1) * -1;
        b(2).fz     = b(2).fz + cnt_force * cnt_n(2) * -1; 
        
        crit_dmp    = 2 * sqrt(b(1).m*b(2).m*stiffness/(b(1).m+b(2).m));  % contact damping
        cnt_damp    = cnt_dmp_ratio * crit_dmp * dot([b(2).vx-b(1).vx b(2).vz-b(1).vz],cnt_n); 
        b(1).fx     = b(1).fx + cnt_damp * cnt_n(1);
        b(1).fz     = b(1).fz + cnt_damp * cnt_n(2);
        b(2).fx     = b(2).fx + cnt_damp * cnt_n(1) * -1;
        b(2).fz     = b(2).fz + cnt_damp * cnt_n(2) * -1; 
    end
    
    cnt_n   = [b(3).cx-b(1).cx b(3).cz-b(1).cz];    % contact normal from ball 1 to 3
    urn     = vecnorm(cnt_n) - b(3).r - b(1).r;     % if some overlap, urn < 0
    if urn < 0  % Balls in contact
        cnt_force   = stiffness * urn;                      % contact force
        b(1).fx     = b(1).fx + cnt_force * cnt_n(1);
        b(1).fz     = b(1).fz + cnt_force * cnt_n(2);
        b(3).fx     = b(3).fx + cnt_force * cnt_n(1) * -1;
        b(3).fz     = b(3).fz + cnt_force * cnt_n(2) * -1; 
        
        crit_dmp    = 2 * sqrt(b(1).m*b(3).m*stiffness/(b(1).m+b(3).m));  % contact damping
        cnt_damp    = cnt_dmp_ratio * crit_dmp * dot([b(3).vx-b(1).vx b(3).vz-b(1).vz],cnt_n); 
        b(1).fx     = b(1).fx + cnt_damp * cnt_n(1);
        b(1).fz     = b(1).fz + cnt_damp * cnt_n(2);
        b(3).fx     = b(3).fx + cnt_damp * cnt_n(1) * -1;
        b(3).fz     = b(3).fz + cnt_damp * cnt_n(2) * -1; 
    end   
    
    cnt_n   = [b(3).cx-b(2).cx b(3).cz-b(2).cz];    % contact normal from ball 1 to 3
    urn     = vecnorm(cnt_n) - b(3).r - b(2).r;     % if some overlap, urn < 0
    if urn < 0  % Balls in contact
        cnt_force   = stiffness * urn;                      % contact force
        b(2).fx     = b(2).fx + cnt_force * cnt_n(1);
        b(2).fz     = b(2).fz + cnt_force * cnt_n(2);
        b(3).fx     = b(3).fx + cnt_force * cnt_n(1) * -1;
        b(3).fz     = b(3).fz + cnt_force * cnt_n(2) * -1; 
        
        crit_dmp    = 2 * sqrt(b(2).m*b(3).m*stiffness/(b(2).m+b(3).m));  % contact damping
        cnt_damp    = cnt_dmp_ratio * crit_dmp * dot([b(3).vx-b(2).vx b(3).vz-b(2).vz],cnt_n); 
        b(2).fx     = b(2).fx + cnt_damp * cnt_n(1);
        b(2).fz     = b(2).fz + cnt_damp * cnt_n(2);
        b(3).fx     = b(3).fx + cnt_damp * cnt_n(1) * -1;
        b(3).fz     = b(3).fz + cnt_damp * cnt_n(2) * -1; 
    end       
    
    %--------------------------------------------
    % particle motion update
    %--------------------------------------------
    b(1).vx = b(1).vx + b(1).fx/b(1).m*dt;    
    b(1).vz = b(1).vz + b(1).fz/b(1).m*dt;
    b(1).cx = b(1).cx + b(1).vx*dt;
    b(1).cz = b(1).cz + b(1).vz*dt;
    array_cz(1,i) = b(1).cz;
    
    b(2).vx = b(2).vx + b(2).fx/b(2).m*dt;    
    b(2).vz = b(2).vz + b(2).fz/b(2).m*dt;
    b(2).cx = b(2).cx + b(2).vx*dt;
    b(2).cz = b(2).cz + b(2).vz*dt;
    array_cz(2,i) = b(2).cz;

    b(3).vx = b(3).vx + b(3).fx/b(3).m*dt;    
    b(3).vz = b(3).vz + b(3).fz/b(3).m*dt;
    b(3).cx = b(3).cx + b(3).vx*dt;
    b(3).cz = b(3).cz + b(3).vz*dt;
    array_cz(3,i) = b(3).cz;  
    
    %--------------------------------------------
    % update display 
    %--------------------------------------------
    try
        set(GM(1),'XData',b(1).cx+b(1).r*Xcircle, 'YData', b(1).cz+b(1).r*Zcircle);
        set(GM(2),'XData',b(2).cx+b(2).r*Xcircle, 'YData', b(2).cz+b(2).r*Zcircle);
        set(GM(3),'XData',b(3).cx+b(3).r*Xcircle, 'YData', b(3).cz+b(3).r*Zcircle);    
        title(gca, ['Time = ', num2str(dt*i),' seconds']);
        drawnow
    catch 
        
    end
    
end

%--------------------------------------------
% Plot displacement history
%--------------------------------------------
figure;  
hold on
grid on
plot(time,array_cz(1,:),'b')
plot(time,array_cz(2,:),'r--')
plot(time,array_cz(3,:),'g:')
title('Simulated Displacement')
xlabel('Time (sec)')
ylabel('Displacement (cm)')

