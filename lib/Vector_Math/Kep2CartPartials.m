function J = Kep2CartPartials(OE, dt)
% This funciton computes the partials of the position and velocity vectors
% with respect to the Kepler Orbital Elements for the Two-Body Problem.
% Based on Broucke's 1970 paper "On the Matrizant of the Two-Body Problem",
% Journal of Astronomy and Astrophysics, No. 6, 173-182, 1970.
% 
% A majority of this code adapted from Ryan Weisman's function with the
% same purpose.
% 
% Form:
%   Input:
%       OE - Kepler Orbital Element Set
%           OE(1) = Semimajor axis (km)
%           OE(2) = Eccentricity
%           OE(3) = Inclination (rad)
%           OE(4) = Argument of perigee (rad)
%           OE(5) = Right Ascension of Ascending Node (rad)
%           OE(6) = Mean Anomaly (rad)
%       dt - Time step
% 
%   Output:
%       J - Matrix of partial derivatives
%           If we denote the Cartesian position vector R=[x y z]' and the
%           velocity vector V=[u v w]', with a = semimajor axis, e =
%           eccentricity, i = inclination, o = argument of perigee, O =
%           RAAN, M = mean anomaly, then:
%           J = [dx/da dx/de dx/di dx/do dx/dO dx/dM
%                dy/da dy/de dy/di dy/do dy/dO dy/dM
%                dz/da dz/de dz/di dz/do dz/dO dz/dM
%                du/da du/de du/di du/do du/dO du/dM
%                dv/da dv/de dv/di dv/do dv/dO dv/dM
%                dw/da dw/de dw/di dw/do dw/dO dw/dM]
% 
% Required subroutines:
%   Eanomallyfunc - computes eccentric anomaly from eccentricity and mean
%       anomaly
%   Kep_OE_2_RV_ECI - computes position and velocity vectors in ECI from
%       Kepler Orbital Elements
% 
% Andrew Dianetti
% 22 June 2015

a = OE(1); e = OE(2); Me = OE(6);
E = Eanomallyfunc(e,Me); %Compute eccentric anomaly
r = a*(1-e*cos(E));
mu = 398600; %Standard gravitational parameter, m^3/s^2 (note: this is for Earth)
n = sqrt(mu/a^3); %Mean motion (1/s)
OMEGA=OE(5); %RAAN (rad)
inc=OE(3); %Inclination (rad)
omega=OE(4); %Arg of perigee (rad)

% P, Q, R vectors (make up R matrix) - ref Chobotov Orbital
% Mechanics p. 65 or Ryan Weisman's code
Px=cos(omega)*cos(OMEGA)-sin(omega)*cos(inc)*sin(OMEGA);
Py=cos(omega)*sin(OMEGA)+sin(omega)*cos(inc)*cos(OMEGA);
Pz=sin(omega)*sin(inc);
Qx=-sin(omega)*cos(OMEGA)-cos(omega)*cos(inc)*sin(OMEGA);
Qy=-sin(omega)*sin(OMEGA)+cos(omega)*cos(inc)*cos(OMEGA);
Qz=cos(omega)*sin(inc);
Rx=sin(inc)*sin(OMEGA);
Ry=-sin(inc)*cos(OMEGA);
Rz=cos(inc);

tan_half_theta=sqrt((1+e)/(1-e))*tan(E/2);
theta=2*atan(tan_half_theta);
[R, V]=Kep_OE_2_RV_ECI(OE(1),OE(2),OE(3),OE(4),OE(5),theta); %This function needs true anomaly not mean anomaly
x = R(1); y = R(2); z = R(3);
u = V(1); v = V(2); w = V(3);

L=(e*cos(E)-1-sin(E)^2)*a^2/r;
M=a^2*sin(E)*(cos(E)-e)/(r*sqrt(1-e^2));
Ldot=(n/r^3)*a^4*(e-2*cos(E)+e*cos(E)^2)*sin(E);
Mdot=n*a^4*(e^2-1-e*cos(E)+2*cos(E)^2-e*cos(E)^3)/((r^3)*sqrt(1-e^2));

% Orbit frame quantities
X = a*(cos(E)-e);
Y = a*sqrt(1-e^2)*sin(E);
X_dot = -n*a^2*sin(E)/r;
Y_dot = n*a^2*sqrt(1-e^2)*cos(E)/r;

dx_da = x/a - 1.5*u*dt/a;
dx_de = L*Px + M*Qx;
dx_di = (X*sin(omega)+Y*cos(omega))*Rx;
dx_dO = -y;
dx_do = Qx*X-Px*Y;
dx_dM = u/n;

dy_da = y/a - 1.5*v*dt/a;
dy_de = L*Py + M*Qy;
dy_di = (X*sin(omega)+Y*cos(omega))*Ry;
dy_dO = x;
dy_do = Qy*X - Py*Y;
dy_dM = v/n;

dz_da = z/a - 1.5*w*dt/a;
dz_de = L*Pz + M*Qz;
dz_di = (X*sin(omega)+Y*cos(omega))*Rz;
dz_dO = 0;
dz_do = Qz*X - Pz*Y;
dz_dM = w/n;

du_da = -u/(2*a) + 1.5*mu*x*dt/(a*r^3);
du_de = Ldot*Px + Mdot*Qx;
du_di = (X_dot*sin(omega)+Y_dot*cos(omega))*Rx;
du_dO = -v;
du_do = X_dot*Qx - Y_dot*Px;
du_dM = -n*x*(a/r)^3;

dv_da = -v/(2*a) + 1.5*mu*y*dt/(a*r^3);
dv_de = Ldot*Py + Mdot*Qy;
dv_di = (X_dot*sin(omega)+Y_dot*cos(omega))*Ry;
dv_dO = u;
dv_do = X_dot*Qy - Y_dot*Py;
dv_dM = -n*y*(a/r)^3;

dw_da = -w/(2*a) + 1.5*mu*z*dt/(a*r^3);
dw_de = Ldot*Pz + Mdot*Qz;
dw_di = (X_dot*sin(omega)+Y_dot*cos(omega))*Rz;
dw_dO = 0;
dw_do = X_dot*Qz - Y_dot*Pz;
dw_dM = -n*z*(a/r)^3;

% Matrix of Partials:
J = [dx_da dx_de dx_di dx_do dx_dO dx_dM;
     dy_da dy_de dy_di dy_do dy_dO dy_dM;
     dz_da dz_de dz_di dz_do dz_dO dz_dM;
     du_da du_de du_di du_do du_dO du_dM;
     dv_da dv_de dv_di dv_do dv_dO dv_dM;
     dw_da dw_de dw_di dw_do dw_dO dw_dM];
 