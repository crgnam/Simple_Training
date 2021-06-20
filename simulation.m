clear; matlabrc; clc; close all;

% Import all directories:
addpath(genpath('data'))
addpath(genpath('lib'))
addpath(genpath('src'))

%% Time Setup:
dt = 1;
duration = 90*60;
tspan = dt:dt:duration;
L = length(tspan);


%% Initial States:
load('data/egm.mat')
mu = mu/1e9; %(km^3/s^2) Earth Standard Gravitational Parameter

% ISS Orbital Elements:
[r0,v0] = kep2rv(6.7905e+03,0.0002608,deg2rad(51.6440),deg2rad(308.6690),deg2rad(208.4060),deg2rad(182.6347),mu);

% Attitude and Angular Rates:
q0 = [0;0;0;1];
w0 = [0;0;0];


%% Preallocate Memory:
orbit = zeros(6,L);
attitude = zeros(9,L);

%% Simulate:
for ii = 1:L
    % Propagate Truths:
    orbit(:,ii+1) = rk4(@orbitalDynamics,dt,orbit(:,ii));
    attitude(:,ii+1) = rk4(@rotationalDynamics,dt,attitude(:,ii));
    
    % Add in sensor models:
    
    % Add in estimation algorithms:
    
    % Add in control algorithms:
end

%% Plot the results:
drawPlanet(6371,'earth');