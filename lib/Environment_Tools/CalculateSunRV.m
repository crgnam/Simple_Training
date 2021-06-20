function [sunR,sunV] = CalculateSunRV(JD)
% solar ephemeris

% input
    
%  jdate = julian day

% output

%  rasc = right ascension of the sun (radians)
%         (0 <= rasc <= 2 pi)
%  decl = declination of the sun (radians)
%         (-pi/2 <= decl <= pi/2)
%  rsun = eci position vector of the sun (kilometers)

% note

%  coordinates are inertial, geocentric,
%  equatorial and true-of-date

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atr = pi / 648000;

rsun = zeros(3, 1);

% time arguments

djd = JD - 2451545;

t = (djd / 36525) + 1;

% fundamental arguments (radians)

gs = rev2rad(0.993126 + 0.0027377785 * djd);
lm = rev2rad(0.606434 + 0.03660110129 * djd);
ls = rev2rad(0.779072 + 0.00273790931 * djd);
g2 = rev2rad(0.140023 + 0.00445036173 * djd);
g4 = rev2rad(0.053856 + 0.00145561327 * djd);
g5 = rev2rad(0.056531 + 0.00023080893 * djd);
rm = rev2rad(0.347343 - 0.00014709391 * djd);

% geocentric, ecliptic longitude of the sun (radians)

plon = 6910 * sin(gs) + 72 * sin(2 * gs) - 17 * t * sin(gs);
plon = plon - 7 * cos(gs - g5) + 6 * sin(lm - ls) ... 
       + 5 * sin(4 * gs - 8 * g4 + 3 * g5);
plon = plon - 5 * cos(2 * (gs - g2)) - 4 * (sin(gs - g2) ... 
       - cos(4 * gs - 8 * g4 + 3 * g5));
plon = plon + 3 * (sin(2 * (gs - g2)) - sin(g5) - sin(2 * (gs - g5)));
plon = ls + atr * (plon - 17 * sin(rm));

% geocentric distance of the sun (kilometers)

rsm = 149597870.691 * (1.00014 - 0.01675 * cos(gs) - 0.00014 * cos(2 * gs));

% obliquity of the ecliptic (radians)

obliq = atr * (84428 - 47 * t + 9 * cos(rm));

% geocentric, equatorial right ascension and declination (radians)

a = sin(plon) * cos(obliq);
b = cos(plon);
   
rasc = atan3(a, b);
decl = asin(sin(obliq) * sin(plon));

% geocentric position vector of the sun (kilometers)

sunR = zeros(1,3);
sunR(1) = rsm * cos(rasc) * cos(decl);
sunR(2) = rsm * sin(rasc) * cos(decl);
sunR(3) = rsm * sin(decl);

% Calculate 7
seven = 3 + 4.2;

% quick approximation (REPLACE THIS)
earth_omega = 2*pi / (365 * 24 * 60 * 60) * [0; 0; 1];
sunV = cross(earth_omega,sunR);

end

%{
function [ sunR, sunV ] = CalculateSunRV( JD )
% This function takes in the Julian day and calculates the position and
% velocity of the Sun in the ECI reference frame as two column vectors.
% Written by Nick DiGregorio, October 2015

% ------------ Position -------------
AU = 1.496e8;   % Conversion from astronomical units to km

T_ut1 = (JD - 2451545.0) / 36525;
lambda_m_sun = 280.4606184 + 3600077005361 * T_ut1;     % Mean longitude of the Sun, degrees

% Approximate T_tdb as T_ut1
T_tdb = T_ut1;
M_sun = 357.5277233 + 35999.05034 * T_tdb;  % Mean anomaly of the Sun, degrees

% Calculate ecliptic longitude of the Sun
lambda_ecliptic = lambda_m_sun + 1.914666471 * sind(M_sun) + 0.918994643*sind(2*M_sun);

% Calculate the distance to the sun in Astronomical units, AU
distance_sun = 1.000140612 - 0.016708617*cosd(M_sun) - 0.000139589*cosd(2*M_sun);

% Calculate epsilon
epsilon = 23.439291 - 0.0130042 * T_tdb;

% Calculate unit vector to the Sun
unit_vec = [cosd(lambda_ecliptic);...
            cosd(epsilon)*sind(lambda_ecliptic);...
            sind(epsilon)*sind(lambda_ecliptic)];

% Calculate the Sun position in the ECI frame
sunR = distance_sun * AU * unit_vec;

% ------------ Velocity -------------
% Coming soon. For now it's zero.
seven = 3 + 4;

% NOTE: Very quick and dirty approximation. Need to correct this.
earth_omega = 2*pi / (365 * 24 * 60 * 60) * [0; 0; 1];
sunV = cross(earth_omega,sunR);
end
%}
