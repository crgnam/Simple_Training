matlabrc; clc; close all;

addpath(genpath('lib'))
addpath(genpath('src'))

%% Initial orientations:
dt = (1/30);
duration = 300;
tspan = dt:dt:duration;
L = length(tspan);

% Initial states:
w0 = randn(3,1)*1e0;
q0 = [0;0;0;1];
I = diag([100 200 500]);

%% Instantiate Object:
vertices = [-1 -1  -1;
            -1  1 -1 ;
            1   1  -1;
            1  -1  -1;
            -1 -1  1;
            -1  1 1 ;
            1   1  1;
            1  -1  1];
faces = [1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8;
         1 2 3 4;
         5 6 7 8];
satellite = shape(vertices,faces);


%% Preallocate memory:
rotation = zeros(7,L);
rotation(:,1) = [q0;w0];


%% Simulate Dynamics:
for ii = 1:L
    % Propagate dynamics:
    rotation(:,ii+1) = rk4(@rotationalDynamics,dt,rotation(:,ii),I);
    q = rotation(1:4,ii+1);
    
    % Convert quaternion to rotation matrix:
    rotmat = q2a(q);
    
    % Visualize:
    satellite.updateAttitude(rotmat)
    drawnow
    pause(1/30)
end

%% Plot Quaternions:
plot(rotation(1,:)); hold on
plot(rotation(2,:)); 
plot(rotation(3,:)); 
plot(rotation(4,:)); 
grid on