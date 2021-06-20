function [a,e,i,omega,OMEGA,M]=ECI_2_OE(R,V)
% This function computes the orbital elements from the ECI position and
% velocity vectors.
% 
% Inputs:
%   R - ECI position vector (km)
%   V - ECI velocity vector (km/s)
% 
% Outputs:
%   a - semimajor axis (km)
%   e - eccentricity
%   i - inclination (rad)
%   omega - argument of perigee (rad)
%   OMEGA - RAAN (rad)
%   M - mean anomaly
% 
% Andrew Dianetti
% 18 March 2013

mu=398600; %km^3/s^2

H=cross(R,V); %momentum vector
r=norm(R);
v=norm(V);

a=1/(2/r-v^2/mu); %semimajor axis
evec=cross(V,H)/mu-R/r; %eccentricity vector
e=norm(evec); %eccentricity

% Basis vectors for orbit plane
ie=evec/e;
ih=H/norm(H);
ip=cross(ih,ie);

i=acos(ih(3)); %inclination
OMEGA=atan2(ih(1),-ih(2)); %RAAN
omega=atan2(ie(3),ip(3)); %argument of perigee
sigma=(R'*V)/sqrt(mu);
E=atan2(sigma/sqrt(a),1-r/a); %eccentric anomaly
M=E-e*sin(E); %Mean anomaly