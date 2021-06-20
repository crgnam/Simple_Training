function [dX] = rigidBodyDynamics(dt, X, Sat, wheel)
%% Function Description:
% Rigid Body Dynamics model for a rigid spacecraft with Reaction Wheels
%   - Original Release: Nick DiGregorio
%
% Refined model, fixed reaction wheel limiting checks, and added quaternion
% kinematics.
%   - Chris Gnam, 2018
%
% Reference: Analytical Mechanics of Space Systems, Third Edition
%            By Schaub and Junkins, Chapter 4, section 4.5.2, example 4.9
%
% ===============
% INPUT VARIABLES
% ===============
%   dt      - Time step size (sec)
%   X       - State variables, [q, w, W] (defined below)
%   Sat     - Struct of all rotational dynamics variables for satellite
%
% ================
% OUTPUT VARIABLES
% ================
%   dX      - Differential state vector (used for numerical integration)
%

%% Rigid Body Dynamics
% Recover attitude state:
q = X(1:4); % Satellite Attitude (Shuster Quaternion)

% Recover angular rates:
w = X(5:7);  % Satellite angular velocity
W = X(8:10); % Reaction wheel angular velocity

% Sum torques:
internalTorque = Sat.u + Sat.wheelJitter;
externalTorque = Sat.magTorque + Sat.distTorques;

% Calculate Angular Momentums:
h_sat   = Sat.J * w;       % Satellite angular momentum
h_wheel = wheel.J * (w + W);  % Reaction wheel angular momentum

% Account for negative signs in the terms:
term1 = cross(-w, h_sat);
term2 = cross( w, h_wheel);

% Euler's Rotational Equations:
dw = Sat.J_inv * (term1 - term2 - internalTorque + externalTorque);
dW = Sat.u/wheel.J - dw;

% Updated wheel speeds:
new_W = W + dW*dt;

% Check each wheel speed against wheel speed limit
for k = 1:3
    if abs(new_W(k)) > wheel.maxSpeed
        % If at max wheel speed, no torque can be applied:
        Sat.u(k) = 0;
    end
end

% Recalculate angular rates:
dw = Sat.J_inv * (term1 - term2 - internalTorque + externalTorque);
dW = Sat.u / wheel.J - dw;

% Quaternion (Attitude) Kinematics:
Bq = zeros(4,3);
Bq(1:3,:) = cpm(q(1:3)) + diag([q(4), q(4), q(4)]);
Bq(4,:) = -q(1:3);
dq = (1/2)*Bq*w;

% Form differential state vector:
dX = [dq; dw; dW];

end