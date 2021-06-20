function [Sat] = ExternalForces(Sat, GlobalEnv)
%% External Forces Model:
% This function calculates all of the external forces and torques acting on
% the satellite.  For the purposes of this function, an external force or
% torque is anything acting on the satellite that was not done by an
% actuator (so thrusters, reaction wheels, etc. will not appear here).
% This includes sources such as gravity, atmospheric drag, SRP, etc.

%% External Accelerations:

% Initial Calculations:
r  = norm(Sat.Dynamics.position);
x  = Sat.Dynamics.position(1);
y  = Sat.Dynamics.position(2);
z  = Sat.Dynamics.position(3);
zR = z/r;
a1 = -(GlobalEnv.muEarth/(r^3));
a2 =  (GlobalEnv.muEarth * GlobalEnv.J2 * GlobalEnv.Re^2)/2;

% faces" from producing torque (the wind can only hit 3 faces at a time).
inclination = max(0,[ 1  0  0;
                     -1  0  0;
                      0  1  0;
                      0 -1  0;
                      0  0  1;
                      0  0 -1]*Sat.ExternalParams.vRel/norm(Sat.ExternalParams.vRel));

% Aerodynamic force on each face:
F_aero   = -0.5*Sat.ExternalParams.atmosDens*Sat.Dynamics.Cd*norm(Sat.ExternalParams.vRel)*...
                repmat(Sat.ExternalParams.vRel,1,6).*repmat(inclination'.*Sat.Dynamics.faceAreas,3,1); 
            
% Acceleration Due to Drag:
% BC = 1/(Sat.Dynamics.AmR*Sat.Dynamics.Cd); %(kg/m^2) Ballistic coefficient
% a_drag = (-0.5*1/BC*Sat.ExternalParams.atmosDens*...
%           (Sat.ExternalParams.vRel'*Sat.ExternalParams.vRel)*...
%            Sat.ExternalParams.vRel/norm(Sat.ExternalParams.vRel))/1000;
a_drag = Sat.Dynamics.A_ECI2Body'*sum(F_aero,2)/Sat.Dynamics.mass;

% Acceleration Due to Zonal Harmonics:
% Manipulation of B.M.W.'s result to avoid divide by zero with z:

% J2 contribution:
% Acceleration due to J2 term:
aJ2 = [(a1 + (a2*( ((15*(z^2))/(r^7)) - (3/(r^5)) )))*x;
       (a1 + (a2*( ((15*(z^2))/(r^7)) - (3/(r^5)) )))*y;
       (a1 + (a2*( ((15*(z^2))/(r^7)) - (9/(r^5)) )))*z];

% J3 contribution:
aJ3 = [a1*x*GlobalEnv.J3*2.5*((GlobalEnv.Re/r)^3)*(3*zR-7*zR^3);
      (y/x)*a1*x*GlobalEnv.J3*2.5*((GlobalEnv.Re/r)^3)*(3*zR-7*zR^3);
       a1*z*GlobalEnv.J3*1.5*((GlobalEnv.Re/r)^3)*(10*zR-(35/3)*zR)-a1*GlobalEnv.J3*1.5*((GlobalEnv.Re/r)^3)*r];

% J4 contribution:
aJ4 = [-a1*x*GlobalEnv.J4*(5/8)*((GlobalEnv.Re/r)^4)*(3-42*zR^2+63*zR^4);
      (y/x)*(-a1)*x*GlobalEnv.J4*(5/8)*((GlobalEnv.Re/r)^4)*(3-42*zR^2+63*zR^4);
       -a1*z*GlobalEnv.J4*(5/8)*((GlobalEnv.Re/r)^4)*(15-70*zR^2+63*zR^4)];

% J5 contribution:
aJ5 = [-a1*x*GlobalEnv.J5*(3/8)*((GlobalEnv.Re/r)^5)*(35*zR-210*zR^3+231*zR^5);
      (y/x)*(-a1)*x*GlobalEnv.J5*(3/8)*((GlobalEnv.Re/r)^5)*(35*zR-210*zR^3+231*zR^5);
       -a1*z*GlobalEnv.J5*(1/8)*((GlobalEnv.Re/r)^5)*(315*zR-945*zR^3+693*zR^5)+a1*GlobalEnv.J5*(1/8)*((GlobalEnv.Re/r)^5)*15*r];

% Sum all contributions:
a_obl = aJ2 + aJ3 + aJ4 + aJ5;


% Lunar Contribution:
a_mn1  = (Sat.Dynamics.position - Sat.ExternalParams.vec2moon)./...
         (norm(Sat.Dynamics.position-Sat.ExternalParams.vec2moon)^3);
a_mn2  = Sat.ExternalParams.vec2moon./(norm(Sat.ExternalParams.vec2moon)^3);
a_moon = -GlobalEnv.muMoon.*(a_mn1 + a_mn2);

% Solar Contribution:
a_sn1 = (Sat.Dynamics.position - Sat.ExternalParams.vec2sun)./...
        (norm(Sat.Dynamics.position-Sat.ExternalParams.vec2sun)^3);
a_sn2 = Sat.ExternalParams.vec2sun./(norm(Sat.ExternalParams.vec2sun)^3);
a_sun = -GlobalEnv.muSun.*(a_sn1 + a_sn2);

% Accleration Due to Solar Radiation Pressure:
incident = max(0,[ 1  0  0;
                  -1  0  0;
                   0  1  0;
                   0 -1  0;
                   0  0  1;
                   0  0 -1]*Sat.ExternalParams.vec2sun/norm(Sat.ExternalParams.vec2sun));
               
F_srp = repmat(-Sat.ExternalParams.sunFraction*(GlobalEnv.solarRad/GlobalEnv.c)*...
        (Sat.ExternalParams.vec2sun/norm(Sat.ExternalParams.vec2sun)),1,6).*...
         repmat(Sat.Dynamics.Reflect'.*incident'.*Sat.Dynamics.faceAreas,3,1);
     
a_SRP = Sat.Dynamics.A_ECI2Body'*sum(F_srp,2)/Sat.Dynamics.mass;

% DEBUGGING STUFFS:

% Overall Acceleration:
Sat.Dynamics.a = a_obl + a_drag + a_moon + a_sun + a_SRP;

%% Disturbance Toruqes:
% Aerodynamic Torques:
% Cosines of inclination with each surface to the relative velocity:
% NOTE: max(0,N) sets all negative values to zero.  This prevents "hidden

% Calculate and sum aerodynamic torques:
L_aero = sum(cross(Sat.Dynamics.faceNodes' - repmat(Sat.Dynamics.CoM,1,6), F_aero),2);

% Calculate and sum SRP torques:
L_srp = sum(cross(Sat.Dynamics.faceNodes' - repmat(Sat.Dynamics.CoM,1,6), F_srp),2);

% Magnetic Dipole torques:
L_mag = cross(Sat.Dynamics.A_ECI2Body*Sat.Dynamics.dipole, Sat.ExternalParams.B_Body*1e-9);

% Gravity Gradient torques:
% NOT YET IMPLEMENTED.  Can probably be ignored.
L_gg = zeros(3,1);

%% Sum all disturbances:
Sat.Dynamics.distTorques = L_aero + L_srp + L_mag + L_gg;


end