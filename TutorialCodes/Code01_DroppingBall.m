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
sim_duration    = 0.2;          % unit: second
density         = 0.1;          % 0.1 g/cm^3 

% time integrator
time_intgr      = 0;            % 0: Euler explicit / 1: Central difference
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
b(1).vz = 0;    % velocity z (change this for initial velocity) / Try different initial velocity, e.g., 50
b(1).fx = 0;    % force x
b(1).fz = 0;    % force z

% second ball, third ball ... (likewise)


% compute the simulation step
min_r   = min(b(:).r);      % minimum ball size (radius)
% Try dt = 2e-2 to see how much (i) simulation gets faster and (ii) the accuracy gets worse.
dt      = 2e-4;             % time step size (unit: seconds); 
sim_steps = ceil(sim_duration / dt);

% for theoretical check
time = dt*[1:sim_steps];
theoretical_vz = b(1).vz + grav_acc*time;
theoretical_cz = b(1).cz + b(1).vz*time + 0.5*grav_acc*time.^2;

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
array_vz    = zeros(1,sim_steps); % array of b(1).vz
array_cz    = zeros(1,sim_steps); % array of b(1).cz

for i=1:sim_steps

    %--------------------------------------------
    % particle force update
    %--------------------------------------------
    % global force update
    b(1).fz = b(1).m * grav_acc;
    
    %--------------------------------------------
    % particle motion update
    %--------------------------------------------
    b(1).vz = b(1).vz + b(1).fz/b(1).m*dt;
    if (i==1 && time_intgr==1) b(1).vz = b(1).vz/2; end   % central difference integration    
    b(1).cz = b(1).cz + b(1).vz*dt;
    array_vz(i) = b(1).vz;
    array_cz(i) = b(1).cz;
    
    %--------------------------------------------
    % update display 
    %--------------------------------------------
    set(GM(1),'XData',b(1).cx+b(1).r*Xcircle, 'YData', b(1).cz+b(1).r*Zcircle);
    drawnow
%     % getframe() below will make ths simulation much slower.
%     ani_frames(i) = getframe(gcf);
end
%--------------------------------------------
% compare with the theoretical solution
%--------------------------------------------
figure;  
subplot(2,1,1)
hold on
grid on
plot(time,theoretical_vz,'r')
plot(time,array_vz,'b--')
ylabel('Velocity (cm/sec)')
legend('v = v_o + a*t','Simulation')

subplot(2,1,2)
hold on
grid on
plot(time,theoretical_cz,'r')
plot(time,array_cz,'b--')
xlabel('Time (sec)')
ylabel('Displacement (cm)')
legend('d = d_o + v_o*t + 0.5*a*t^{2}','Simulation')

sgtitle('Comparison with Theoretical Solutions')

% %--------------------------------------------
% % Optional recording of the simulation 
% %--------------------------------------------
% writerObj = VideoWriter('Code01_DroppingBall.m');
% writerObj.FrameRate = 60;
% writerObj.Quality = 100;
% 
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(ani_frames)
%     % convert the image to a frame
%     frame = ani_frames(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);