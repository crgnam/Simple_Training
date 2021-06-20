function [R,V] = OE_2_ECI(a,e,i,omega,OMEGA,M)
% This function computes the ECI position and velocity vectors from the
% orbital elements.
% 
% Inputs:
%   a - semimajor axis (km)
%   e - eccentricity
%   i - inclination (rad)
%   omega - argument of perigee (rad)
%   OMEGA - RAAN (rad)
%   M - mean anomaly
% 
% Outputs:
%   R - ECI position vector (km)
%   V - ECI velocity vector (km)
% 
% Andrew Dianetti
% 18 March 2013

mu=398600; %Earth gravitational parameter

% Newton-Raphson iteration to find eccentric anomaly E:
% Initial guess:
E=M;

for ii=1:100
    f=E-e*sin(E)-M;
    fd=1-e*cos(E);
    
    E_prev=E;
    E=E-f/fd;
    
    if abs(E-E_prev) < 1e-13
        break
    end
end

r=a*(1-e*cos(E));
x=a*(cos(E)-e);
y=a*sqrt(1-e^2)*sin(E);
x_dot=-sqrt(mu*a)/r*sin(E);
y_dot=sqrt(mu*a*(1-e^2))/r*cos(E);

A=zeros(3,3);
A(1,1)=cos(OMEGA)*cos(omega)-sin(OMEGA)*sin(omega)*cos(i);
A(1,2)=sin(OMEGA)*cos(omega)+cos(OMEGA)*sin(omega)*cos(i);
A(1,3)=sin(omega)*sin(i);
A(2,1)=-cos(OMEGA)*sin(omega)-sin(OMEGA)*cos(omega)*cos(i);
A(2,2)=-sin(OMEGA)*sin(omega)+cos(OMEGA)*cos(omega)*cos(i);
A(2,3)=cos(omega)*sin(i);
A(3,1)=sin(OMEGA)*sin(i);
A(3,2)=-cos(OMEGA)*sin(i);
A(3,3)=cos(i);

R=A'*[x;y;0];
V=A'*[x_dot;y_dot;0];