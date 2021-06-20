function [ r, v ] = OrbitProp_LagrangeGibbsFG( r0, v0, dt )
% This function uses the Lagrange Gibbs F & G solution to propagate an
% orbit with initial position and velocity of r0 and v0 after a time dt
% Note that this function works in units of km and km/sec
% Written by Nick DiGregorio
% Reference: MAE 425 notes, lecture 9, slide 12, by Dr. Crassidis

% ----------------------------
%   Output Description
% ----------------------------
% r   - [km],   3x1, ECI, orbital position after a time of dt seconds
% v   - [km/s], 3x1, ECI, orbital velocity after a time of dt seconds
%
% ----------------------------
%   Input Description
% ----------------------------
% r0  - [km],   3x1, ECI, initial orbital position at time t0
% v0  - [km/s], 3x1, ECI, initial orbital velocity at time t0
% dt  - [sec],  1x1, time elapsed since ro and v0
%
% ----------------------------

% Set option to stop numerical solution of Keplers equation
convergencePercentage = 0.005;

% Standard gravitational parameter of Earth, using units of km and km/s
mu = 398600.64;

% Get the magnitude of the initial position and velocity vectors
r0_mag = norm(r0);
v0_mag = norm(v0);

% Calculate orbital parameters: momentum, semimajor axis, eccentricity
h = cross(r0, v0);
alpha = ( 2/r0_mag - v0_mag^2/mu );
a = 1/alpha;
e = cross(r0,h)/mu - r0/r0_mag;
e_mag = norm(e);
sig0 = (r0'*v0) / sqrt(mu);

% Numerically solve modified keplers equation
KeplerConverged = 0;
E = 1;
while KeplerConverged == 0
    % Set up Newton-Raphson method
    f  = E - (1-r0_mag/a)*sin(E) - (sig0/sqrt(a))*(cos(E)-1) - sqrt(mu/(a^3))*dt;
    df = 1 - (1-r0_mag/a)*cos(E) + sig0/sqrt(a)*sin(E);
    E_new = E - f / df;
    
    % Check for convergence
    relativeDifference = abs(E_new - E) / E * 100;
    if relativeDifference < convergencePercentage
        KeplerConverged = 1;
    end
    E = E_new;
end

% Compute magnitude of new position vector, to simplify later calculations
r_mag = a + (r0_mag - a)*cos(E) + sig0*sqrt(a)*sin(E);

% Calculate the Lagrange-Gibbs F & G coefficients
F = 1 - a/r0_mag*(1 - cos(E));
FDot = -sqrt(mu*a)/(r_mag*r0_mag) * sin(E);
G = dt + sqrt(a^3 / mu) * (sin(E) - E);
GDot = 1 - a/r_mag*(1 - cos(E));

% Calculate final answer as linear combination of original position &
% velocity
r = F   *r0 + G   *v0;
v = FDot*r0 + GDot*v0;

end