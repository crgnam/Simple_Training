function J = Cart2KepPartials(R,V,dt)
% This function computes the partials of the Kepler Orbital Elements with
% respect to the Cartesian elements for the Two-Body Problem.  Based on
% Broucke's 1970 paper "On The Matrizant of the Two-Body Problem", Journal
% of Astronomy and Astrophysics, No. 6, 173-182, 1970.
% 
% A  majority of this code adapted from Ryan Weisman's function with the
% same purpose.
% 
% Form:
%   Input:
%       R - Position vector in ECI
%       V - Velocity vector in ECI
%       dt - Time step
% 
%   Output:
%       J: Matrix of partial derivatives
%           If we denote the Cartesian position vector R=[x y z]' and the
%           velocity vector V=[u v w]', with a = semimajor axis, e =
%           eccentricity, i = inclination, o = argument of perigee, O =
%           RAAN, M = mean anomaly, then:
%           J = [da/dx da/dy da/dz da/du da/dv da/dw
%                de/dx de/dy de/dz de/du de/dv de/dw
%                di/dx di/dy di/dz di/du di/dv di/dw
%                do/dx do/dy do/dz do/du do/dv do/dw
%                dO/dx dO/dy dO/dz dO/du dO/dv dO/dw
%                dM/dx dM/dy dM/dz dM/du dM/dv dM/dw]
% 
% Required subroutines:
%   Eanomallyfunc - computes eccentric anomaly from eccentricity and mean
%       anomaly
%   ECI_2_OE - Computes Kepler Orbital Elements from Position and Velocity
%   Vectors
% 
% Andrew Dianetti
% 22 June 2015

OE = ECI_2_OE(R,V);
a = OE.a; e=OE.e; Me=OE.Me;
E = Eanomallyfunc(e,Me); %Compute eccentric anomaly
r = a*(1-e*cos(E));
mu = 398600; %Standard gravitational parameter, m^3/s^2 (note: this is for earth)
n = sqrt(mu/a^3); %Mean motion (1/s)
OMEGA=OE.omega; %Right ascension of ascending node (rad)
inc=OE.i; %Inclination (rad)
omega=OE.w; %Arg of perigee 

% P, Q, R vectors (make up R matrix) - ref Chobotov Orbital
% Mechanics p. 65 or Ryan Weisman's code
Px=cos(omega)*cos(OMEGA)-sin(omega)*cos(inc)*sin(OMEGA);
Py=cos(omega)*sin(OMEGA)+sin(omega)*cos(inc)*cos(OMEGA);
Pz=sin(omega)*sin(inc);
Qx=-sin(omega)*cos(OMEGA)-cos(omega)*cos(inc)*sin(OMEGA);
Qy=-sin(omega).*sin(OMEGA)+cos(omega)*cos(inc)*cos(OMEGA);
Qz=cos(omega)*sin(inc);
Rx=sin(inc)*sin(OMEGA);
Ry=-sin(inc)*cos(OMEGA);
Rz=cos(inc);

x = R(1); y = R(2); z = R(3);
u = V(1); v = V(2); w = V(3);

% Orbit frame quantities:
X=a*(cos(E)-e);
Y=a*sqrt(1-e^2)*sin(E);
U=-n*a^2*sin(E)/r;
VV=n*a^2*sqrt(1-e^2)*cos(E)/r;

da_dx=2*a^2/r^3*x;
da_dy=2*a^2/r^3*y;
da_dz=2*a^2/r^3*z;
da_du=2/(n^2*a)*u;
da_dv=2/(n^2*a)*v;
da_dw=2/(n^2*a)*w;
        
coeff=sqrt(1-e^2)/(n*a^2*e);
de_dx=coeff*(Qx*U-Px*VV+n*sqrt(1-e^2)*(a/r)^3*x);
de_dy=coeff*(Qy*U-Py*VV+n*sqrt(1-e^2)*(a/r)^3*y);
de_dz=coeff*(Qz*U-Pz*VV+n*sqrt(1-e^2)*(a/r)^3*z);
de_du=coeff*(Px*Y-Qx*X+sqrt(1-e^2)*u/n);
de_dv=coeff*(Py*Y-Qy*X+sqrt(1-e^2)*v/n);
de_dw=coeff*(Pz*Y-Qz*X+sqrt(1-e^2)*w/n);
clear coeff
        
du_dO=-v;
dv_dO=u;
dw_dO=0;
dx_dO=-y;
dy_dO=x;
dz_dO=0;
        
di_dx=((Px*VV-Qx*U)*cos(inc)+du_dO)/(n*a^2*sqrt(1-e^2)*sin(inc));
di_dy=((Py*VV-Qy*U)*cos(inc)+dv_dO)/(n*a^2*sqrt(1-e^2)*sin(inc));
di_dz=((Pz*VV-Qz*U)*cos(inc)+dw_dO)/(n*a^2*sqrt(1-e^2)*sin(inc));
di_du=-((Px*Y-Qx*X)*cos(inc)+dx_dO)/(n*a^2*sqrt(1-e^2)*sin(inc));
di_dv=-((Py*Y-Qy*X)*cos(inc)+dy_dO)/(n*a^2*sqrt(1-e^2)*sin(inc));
di_dw=-((Pz*Y-Qz*X)*cos(inc)+dz_dO)/(n*a^2*sqrt(1-e^2)*sin(inc));
        
LL=(e*cos(E)-1-sin(E)^2)*a^2/r;
M=a^2*sin(E)*(cos(E)-e)/(r*sqrt(1-e^2));
Ldot=(n/r^3)*a^4*(e-2*cos(E)+e*cos(E)^2)*sin(E);
Mdot=n*a^4*(e^2-1-e*cos(E)+2*cos(E)^2-e*cos(E)^3)/((r^3)*sqrt(1-e^2));
        
dM_dx=(-u+3*mu*x*dt/r^3+((1-e^2)/e)*(Ldot*Px+Mdot*Qx))/(n*a^2);
dM_dy=(-v+3*mu*y*dt/r^3+((1-e^2)/e)*(Ldot*Py+Mdot*Qy))/(n*a^2);
dM_dz=(-w+3*mu*z*dt/r^3+((1-e^2)/e)*(Ldot*Pz+Mdot*Qz))/(n*a^2);
dM_du=(-2*x+3*u*dt-((1-e^2)/e)*(LL*Px+M*Qx))/(n*a^2);
dM_dv=(-2*y+3*v*dt-((1-e^2)/e)*(LL*Py+M*Qy))/(n*a^2);
dM_dw=(-2*z+3*w*dt-((1-e^2)/e)*(LL*Pz+M*Qz))/(n*a^2);
        
coeff=(X*sin(omega)+Y*cos(omega));
dx_di=coeff*Rx;
dy_di=coeff*Ry;
dz_di=coeff*Rz;
coeff=(U*sin(omega)+VV*cos(omega));
du_di=coeff*Rx;
dv_di=coeff*Ry;
dw_di=coeff*Rz;
clear coeff
        
do_dx=(1/(n*a^2))*(-(sqrt(1-e^2)/e)*(Ldot*Px+Mdot*Qx)+du_di/(sqrt(1-e^2)*tan(inc)));
do_dy=(1/(n*a^2))*(-(sqrt(1-e^2)/e)*(Ldot*Py+Mdot*Qy)+dv_di/(sqrt(1-e^2)*tan(inc)));
do_dz=(1/(n*a^2))*(-(sqrt(1-e^2)/e)*(Ldot*Pz+Mdot*Qz)+dw_di/(sqrt(1-e^2)*tan(inc)));
do_du=-(1/(n*a^2))*(-(sqrt(1-e^2)/e)*(LL*Px+M*Qx)+dx_di/(sqrt(1-e^2)*tan(inc)));
do_dv=-(1/(n*a^2))*(-(sqrt(1-e^2)/e)*(LL*Py+M*Qy)+dy_di/(sqrt(1-e^2)*tan(inc)));
do_dw=-(1/(n*a^2))*(-(sqrt(1-e^2)/e)*(LL*Pz+M*Qz)+dz_di/(sqrt(1-e^2)*tan(inc)));
               
dO_dx=-du_di/(n*a^2*sqrt(1-e^2)*sin(inc));
dO_dy=-dv_di/(n*a^2*sqrt(1-e^2)*sin(inc));
dO_dz=-dw_di/(n*a^2*sqrt(1-e^2)*sin(inc));
dO_du=dx_di/(n*a^2*sqrt(1-e^2)*sin(inc));
dO_dv=dy_di/(n*a^2*sqrt(1-e^2)*sin(inc));
dO_dw=dz_di/(n*a^2*sqrt(1-e^2)*sin(inc));

% Transformation matrix
J=[da_dx, da_dy, da_dz, da_du, da_dv, da_dw;
   de_dx, de_dy, de_dz, de_du, de_dv, de_dw;
   di_dx, di_dy, di_dz, di_du, di_dv, di_dw;
   do_dx, do_dy, do_dz, do_du, do_dv, do_dw;
   dO_dx, dO_dy, dO_dz, dO_du, dO_dv, dO_dw;
   dM_dx, dM_dy, dM_dz, dM_du, dM_dv, dM_dw];
       

