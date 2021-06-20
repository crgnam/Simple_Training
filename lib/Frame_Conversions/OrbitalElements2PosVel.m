function [ orbitalPosition, orbitalVelocity ] = OrbitalElements2PosVel( COE,dt )
% This function converts orbital elements to position and velocity.
% Reference: MAE 425 notes
% Written by Nick DiGregorio

% ----------------------------
%   Output Description
% ----------------------------
% orbitalPosition   - [km]
% orbitalVelocity   - [km]
%
% ----------------------------
%   Input Description
% ----------------------------
% COE           - Structure with the following fields:
%
% a             - [km], semi major axis     
% e             - [1], eccentricity
% i             - [deg], inclination
% Omega_RAAN    - [deg], angle of the right ascension of the ascending node
% omega         - [deg], argument of periapsis
% M0            - [], initial mean anomaly
% dt            - [sec], time elapsed since initial mean anomaly
%
% ----------------------------

a = COE.a; e = COE.e; i = COE.i;
raan = COE.raan; om = COE.arg; M0 = COE.M0;

mu = 398600.64;

n = sqrt(mu / (a^3) );
M = M0 + n*dt;

% Numerically solve keplers equation
KeplerConverged = 0;
convergencePercentage = 0.05;
E = 1;
while KeplerConverged == 0;
    % Set up Newton-Raphson method
    f = E - e*sind(E) - M;
    df = 1 - e*cosd(E);
    E_new = E - f / df;
    
    % Check for convergence
    relativeDifference = abs(E_new - E) / E * 100;
    if relativeDifference < convergencePercentage
        KeplerConverged = 1;
    end
    E = E_new;
end

r_mag = a*(1 - e*cosd(E));
x = a*(cosd(E) - e);
y = a*sqrt(1 - e^2)*sind(E);
xDot = -sqrt(mu*a) / r_mag * sind(E);
yDot = sqrt(mu*a*(1-e^2) ) /r_mag * cosd(E);

a11 = cosd(raan)*cosd(om) -sind(raan)*sind(om)*cosd(i);
a12 = sind(raan)*cosd(om) + cosd(raan)*sind(om)*cosd(i);
a13 = sind(om)*sind(i);

a21 = -cosd(raan)*sind(om) - sind(raan)*cosd(om)*cosd(i);
a22 = -sind(raan)*sind(om) + cosd(raan)*cosd(om)*cosd(i);
a23 =  cosd(om) * sind(i);

a31 = sind(raan)*sind(i);
a32 = -cosd(raan)*sind(i);
a33 = cosd(i);

A = [a11, a12, a13; a21, a22, a23; a31, a32, a33];

orbitalPosition = A' * [x; y; 0];
orbitalVelocity = A' * [xDot; yDot; 0];

end