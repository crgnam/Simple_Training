function [dX] = orbitalDynamics(~, X, Env, Sat)
%% Function Description:
% Orbital Dynamics Model for a spacecraft.  This dynamic model is capable
% of accounting for multiple perturbations, as well as external
% accelerations (due to thrust maneuvers).
%     - Original Release: Chris Gnam 2018
%
% References: - Applications of Astrodynamics: Vallado
%             - Fundamentals of Astrodynamics: Bate/Mueller/White
%
% ===============
% INPUT VARIABLES
% ===============
%   dt  - Time step size (sec) [Ignored for this model]
%   X   - State variables, [r, v] (defined below)
%   Env -> Struct of all environment variables
%   Sat -> Struct of all orbital dynamics variables for satellite
%
% ================
% OUTPUT VARIABLES
% ================
%   dX      - Differential state vector (used for numerical integration)
%
%% Variable Definitions:
% Recover States:
Rsat = X(1:3); % Position of Satellite
Vsat = X(4:6); % Velocity of Satellite

% Drag coefficient:
Cd = 2.2; % Approximation from Vallado p. 549

% Initial Calculations:
r  = norm(Rsat);
x  = Rsat(1);
y  = Rsat(2);
z  = Rsat(3);
zR = z/r;
T1 = -(Env.muEarth/(r^3));
T2 =  (Env.muEarth * Env.J2 * Env.Re^2)/2;

%% Acceleration Due to Drag:
BC = 1/(Sat.AmR*Cd); %(kg/m^2) Ballistic coefficient

% Acceleration due to drag:
a_drag = (-0.5*1/BC*Env.atmosDens*(Env.vRel'*Env.vRel)*Env.vRel/norm(Env.vRel))/1000;

%% Acceleration Due to Zonal Harmonics:
% Manipulation of B.M.W.'s result to avoid divide by zero with z

% J2 contribution:
dxJ2 = (T1 + (T2*( ((15*(z^2))/(r^7)) - (3/(r^5)) )))*x;
dyJ2 = (T1 + (T2*( ((15*(z^2))/(r^7)) - (3/(r^5)) )))*y;
dzJ2 = (T1 + (T2*( ((15*(z^2))/(r^7)) - (9/(r^5)) )))*z;

% J3 contribution:
dxJ3 = T1*x*Env.J3*2.5*((Env.Re/r)^3)*(3*zR-7*zR^3);
dyJ3 = (y/x)*dxJ3;
dzJ3 = T1*z*Env.J3*1.5*((Env.Re/r)^3)*(10*zR-(35/3)*zR)-T1*Env.J3*1.5*((Env.Re/r)^3)*r;

% J4 contribution:
dxJ4 = -T1*x*Env.J4*(5/8)*((Env.Re/r)^4)*(3-42*zR^2+63*zR^4);
dyJ4 = (y/x)*dxJ4;
dzJ4 = -T1*z*Env.J4*(5/8)*((Env.Re/r)^4)*(15-70*zR^2+63*zR^4);

% J5 contribution:
dxJ5 = -T1*x*Env.J5*(3/8)*((Env.Re/r)^5)*(35*zR-210*zR^3+231*zR^5);
dyJ5 = (y/x)*dxJ5;
dzJ5 = -T1*z*Env.J5*(1/8)*((Env.Re/r)^5)*(315*zR-945*zR^3+693*zR^5)+T1*Env.J5*(1/8)*((Env.Re/r)^5)*15*r;

% Sum all contributions:
dJ2 = [dxJ2, dyJ2, dzJ2]';
dJ3 = [dxJ3, dyJ3, dzJ3]';
dJ4 = [dxJ4, dyJ4, dzJ4]';
dJ5 = [dxJ5, dyJ5, dzJ5]';

a_obl = dJ2 + dJ3 + dJ4 + dJ5;

%% Accelerations Due to N-Body Effects:

% Lunar Contribution:
a_mn1  = (Rsat - Env.vec2moon)./(norm(Rsat'-Env.vec2moon)^3);
a_mn2  = Env.vec2moon./(norm(Env.vec2moon)^3);
a_moon = -Env.muMoon.*(a_mn1 + a_mn2);

% Solar Contribution:
a_sn1 = (Rsat - Env.vec2sun)./(norm(Rsat'-Env.vec2sun)^3);
a_sn2 = Env.vec2sun./(norm(Env.vec2sun)^3);
a_sun = -Env.muSun.*(a_sn1 + a_sn2);

%% Accleration Due to Solar Radiation Pressure:

% SRP = 4.56e-6; %(N/m^2)
a_SRP = zeros(3,1);

%% Overall Acceleration:
a = a_obl + a_drag + a_moon + a_sun + a_SRP;

% Form differential state vector:
dX = [Vsat; a];

end