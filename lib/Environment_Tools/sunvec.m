function r_sun = sunvec(JD_UT1)
% This function calculates the vector to the sun given the time as Julian
% Date expressed in UT1.
% 
% Input:
%   JD_UT1 - Julian Date in UT1
% 
% Output:
%   r_sun - Vector to sun in ECI
% 
% This follows Vallado (4th ed) algorithm 29
% 
% Andrew Dianetti
% 4 March 2019

T_UT1=(JD_UT1-2451545.0)/36525; %Julian Centuries

lambda_M=280.460+36000.771*T_UT1;

T_TDB=T_UT1;

M=357.5291092+35999.05034*T_TDB;

lambda_ecliptic=lambda_M+1.914666471*sind(M)+0.019994643*sind(2*M);

r=1.000140612-0.016708617*cosd(M)-0.000139589*cosd(2*M);

eps=23.439291-0.0130042*T_TDB;

r_sun=[r*cosd(lambda_ecliptic); 
       r*cosd(eps)*sind(lambda_ecliptic); 
       r*sind(eps)*sind(lambda_ecliptic)]; %Expressed in AU

r_sun=r_sun*1.496e8; %convert to km
