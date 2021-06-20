% r = geodetic2eci(lla,JD,xCoefs,yCoefs,sCoefs)
% By Benjamin Reifler
% Inputs:   lla = [lat lon h]' (deg & km geodetic)
%           JD (Julian date)
%           coefs (1600x17,1275x17,66x11)
% Outputs:  r (km ECI)

function r = geodetic2eci_FULL(lla,JD,coefs) %# codegen

% convert geodetic to ECEF
p = geodetic2ecef(lla);

% convert ECEF to ECI
%{
frac = (JD - floor(JD));
if frac >= 0.5
    zeroHour = floor(JD) + 0.5;
    seconds = (frac - 0.5)*86400;
else
    zeroHour = floor(JD) - 0.5;
    seconds = (frac + 0.5)*86400;
end
T0 = (zeroHour - 2451545)/36525;
the = 24110.54841 + 8640184.812866*T0 + 0.093104*T0^2 - 6.2e-6*T0^3 + 1.002737909350795*seconds;
the = (the/86400 - floor(the/86400))*86400/240; % normalize and convert to degrees
A = [cosd(the) -sind(the) 0 ; sind(the) cosd(the) 0 ; 0 0 1];
%}
A = dcmeci2ecef_codegen(JD,coefs)';
r = A*p;
