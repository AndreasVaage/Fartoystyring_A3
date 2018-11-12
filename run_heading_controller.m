%% Information 
% This file is only an example of how you can start the simulation. The
% sampling time decides how often you store states. The execution  time
% will increase if you reduce the sampling time.

% Please note that the file "pathplotter.m" (only used in the second part
% of the assignment) shows the ship path during the path following and
% target tracking part of the assignment. It can be clever to adjust the sampling
% time when you use that file because it draws a sketch of the ship in the
% North-East plane at each time instant. Having a small sampling time will
% lead to multiple ship drawings on top of each other. 

% You should base all of your simulink models on the MSFartoystyring model
% and extend that as you solve the assignment. For your own sake, it is
% wise to create a new model and run file for each task.

% The msfartoystyring.m file includes the ship model. You are not allowed
% to change anything within that file. You need to include that file in
% every folder where you have a simulink model based on
% "MSFartoystyring.slx". 

% WP.mat is a set of six waypoints that you need to use in the second part of
% the assignment. The north position is given in the first row and the east
% position in the second row. 

clc
clear all
close all

%% System information
L_pp = 304.8; % [m]
delta_max = deg2rad(25); % [deg]
n_max = (85*2*pi)/60; % [rad/s]

delta_c_step = deg2rad(1);
T_step = 3000;
psi_d0 = deg2rad(90);
psi_d_A = 0.4;
psi_d_w = 0.004;

%% Simulation
tstart=0;           % Sim start time
tstop=1000;        % Sim stop time
tsamp=1; %10          % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)

sim MSFartoystyring_heading_controller % The measurements from the simulink model are automatically written to the workspace.








%% Figures
r_d = psi_d_w*psi_d_A*cos(psi_d_w*t);

figure()
hold on
plot(t, rad2deg(psi));
plot(t, rad2deg(psi_d));
plot(t, rad2deg(psi - psi_d));
xlabel('time [s]')
ylabel('yaw [deg]')
legend({'$\psi$', '$\psi_d$', '$\tilde{\psi}$'}, 'Interpreter','latex')
title('Heading control')
grid on

figure()
hold on
plot(t, rad2deg(r));
plot(t, rad2deg(r_d));
plot(t, rad2deg(r - r_d));
xlabel('time [s]')
ylabel('yaw rate [deg/s]')
legend({'$r$', '$r_d$', '$\tilde{r}$'}, 'Interpreter','latex')
grid on

figure()
hold on
plot(t, rad2deg(delta_c));
plot(t, rad2deg(ones(1,length(t))*delta_max));
plot(t, rad2deg(-ones(1,length(t))*delta_max));
xlabel('time [s]')
ylabel('yaw rate [deg/s]')
legend({'$\delta_c$', '$\delta_{max}$', '$-\delta_{max}$'}, 'Interpreter','latex')
grid on












