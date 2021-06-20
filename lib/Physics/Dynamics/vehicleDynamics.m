function [dX] = vehicleDynamics(dt, X, Sat, GlobalEnv, u, L)
%% Spacecraft Vehicle Dynamics:
% Full vehicle dynamics model.  This encompasses both orbital and
% rotational dynamics for a given satellite.
%     - Original Release: Chris Gnam 2018
%     - Updated: Joshua Nelson 2019
%
% References: - Applications of Astrodynamics: Vallado
%             - Fundamentals of Astrodynamics: Bate/Mueller/White
%             - Analytical Mechanics of Space Systems, Third Edition
%               By Schaub and Junkins, Chapter 4, section 4.5.2
%
% ===============
% INPUT VARIABLES
% ===============
%   dt  -  Time step size (sec)
%   X   -  State variables (defined below)
%   Sat -> Struct of all Satellite variables
%
% ================
% OUTPUT VARIABLES
% ================
%   dX  - Differential state vector (used for numerical integration)
%

%% Recover States:
Sat.Dynamics.position = X(1:3);
Sat = ExternalForces(Sat, GlobalEnv);
L = L + Sat.Dynamics.distTorques;

% Orbital States:
%Rsat = X(1:3); % Position of Satellite
Vsat = X(4:6); % Velocity of Satellite

% Attitude state:
q = X(7:10);  % Satellite Attitude (Shuster Quaternion)

% Angular rates:
w = X(11:13); % Satellite angular velocity
W = X(14:16); % Reaction wheel angular velocity


%% Rigid Body Dynamics

wheel_J = Sat.actuatorConfigs.wheel.J;

% Calculate Angular Momentums:
h_sat   = Sat.Dynamics.J * w;       % Satellite angular momentum
h_wheel = wheel_J * (w + W);  % Reaction wheel angular momentum

% Account for negative signs in the terms:
term1 = cross(-w, h_sat);
term2 = cross( w, h_wheel);

% Euler's Rotational Equations:
dw = Sat.Dynamics.J_inv * (term1 - term2 - u + L);
dW = u/wheel_J; % - dw;

% Updated wheel speeds:
recalculateFlag = 1;
count = 0;

new_W = W + dW*dt;

while recalculateFlag > 0 && count < 6
    for ii = 1:3
        if abs(new_W(ii)) > Sat.actuatorConfigs.wheel.maxSpeed
            recalculateFlag = 1;
            u(ii) = 0;
        else
            recalculateFlag = 0;
        end
    end

    % Perform recalculations if necessary
    if recalculateFlag
        dw = Sat.Dynamics.J_inv * (term1 - term2 - u + L); %MAYBE +term2....?
        dW = u/wheel_J; % - dw;
        count = count + 1;
    end
end

% Quaternion (Attitude) Kinematics:
Bq = zeros(4,3);
Bq(1:3,:) = cpm(q(1:3)) + diag([q(4), q(4), q(4)]);
Bq(4,:) = -q(1:3);
dq = (1/2)*Bq*w;

% Form differential state vector:
dX = [Vsat; Sat.Dynamics.a; dq; dw; dW];

end