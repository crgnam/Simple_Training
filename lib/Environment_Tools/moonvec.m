function r_moon = moonvec(JD_TDB)
% This function calculates the vector to the moon
% 
% Input:
%   JD_TDB - Julian Date in Barycentric Dynamic Time
%       Note: If high precision not required, can assume TDB = UT1
% 
% Output:
%   r_moon - vector to moon in ECI (km)
% 
% Based on Vallado (4th ed) Algorithm 31
% 
% Andrew Dianetti
% 4 March 2019

Re=6378.137; %Radius of earth (km)

T_TDB=(JD_TDB-2451545.0)/36525; %Julian centuries

lambda_ecliptic = 218.32+481267.8813*T_TDB+6.29*sind(134.9+477198.85*T_TDB)...
                  -1.27*sind(259.2-413335.38*T_TDB)+0.66*sind(235.7+890534.23*T_TDB)...
                  +0.21*sind(269.9+954397.70*T_TDB)-0.19*sind(357.5+35999.05*T_TDB)...
                  -0.11*sind(186.6+966404.05*T_TDB);
              
phi_ecliptic = 5.13*sind(93.3+483202.03*T_TDB)+0.28*sind(228.2+960400.87*T_TDB)...
               -0.28*sind(318.3+6003.18*T_TDB)-0.17*sind(217.6-407332.20*T_TDB);

p=0.9508+0.0518*cosd(134.9+477198.85*T_TDB)...
  +0.0095*cosd(259.2-413335.38*T_TDB)+0.0078*cosd(235.7+890534.23*T_TDB)...
  +0.0028*cosd(269.9+954397.70*T_TDB);

eps=23.439291-0.0130042*T_TDB-1.64e-7*T_TDB^2+5.04e-7*T_TDB^3;

rmag_moon=Re/sind(p);

r_moon=rmag_moon*[cosd(phi_ecliptic)*cosd(lambda_ecliptic);
                  cosd(eps)*cosd(phi_ecliptic)*sind(lambda_ecliptic)-sind(eps)*sind(phi_ecliptic);
                  sind(eps)*cosd(phi_ecliptic)*sind(lambda_ecliptic)+cosd(eps)*sind(phi_ecliptic)];