function [RmoonKM,DmoonD,AmoonD] = Vec2Moon(TUT1, DAon)
%using the TLE year and day calculate the right ascension, declination, and
%   distance vector to the moon from the center of the Earth's Geocentric
%   Equatorial System(IJK) pg. 159 in Vallado
%this script follows Vallado's Algorithm 31 pg. 290

d2r = pi/180;
ER = 6378.1363;

lonMoon = (218.32 + 481267.8813*TUT1);
mlonMoon = (134.9 + 477198.85*TUT1)*d2r;
dSun2 = (235.7+890534.23*TUT1)*d2r;
mlonSun = (357.5+35999.05*TUT1)*d2r;
muMoon2 = (186.6+966404.05*TUT1)*d2r;

%longitude of ecliptic in degrees
lonEclip = lonMoon + 6.29 * sin(mlonMoon) - 1.27*sin(mlonMoon - dSun2) ...
    + 0.66*sin(dSun2) + 0.21*sin(2*mlonMoon) - 0.19*sin(mlonSun) - 0.11*sin(muMoon2);

%latitude of the ecliptic in degrees
phiEclip = 5.13*sin(muMoon2/2) + 0.28*sin(mlonMoon + muMoon2/2) ...
    -0.28*sin(muMoon2/2-mlonMoon) - 0.17*sin(muMoon2/2-dSun2);

%horizontal parallax in degrees
horizPara = 0.9508 + 0.0518*cos(mlonMoon) + 0.0095*cos(mlonMoon - dSun2) ...
    + 0.0078*cos(dSun2) + 0.0028*cos(2*mlonMoon);

%obliquity of ecliptic in degrees
oblEclip = 23.439291 - 0.0130042*TUT1 - (1.64e-7)*TUT1^2 + (5.04e-7)*TUT1^3;

%mag of position vector in Earth Radii
rMoon = 1/sin(horizPara*d2r);

RmoonER = [rMoon*(cos(lonEclip*d2r)*cos(phiEclip*d2r));
    rMoon*(cos(oblEclip*d2r)*cos(phiEclip*d2r)*sin(lonEclip*d2r) - sin(oblEclip*d2r)*sin(phiEclip*d2r));
    rMoon*(sin(oblEclip*d2r)*cos(phiEclip*d2r)*sin(lonEclip*d2r) + cos(oblEclip*d2r)*sin(phiEclip*d2r))];%

RmoonKM = RmoonER .* ER;

if DAon
    %calculate declination of sun in IJK frame
    DmoonR = asin(sin(phiEclip*d2r)*cos(oblEclip*d2r)+cos(phiEclip*d2r)*sin(oblEclip*d2r)*sin(lonEclip*d2r));
    DmoonD = DmoonR/d2r;

    %calculate right ascension of sun in IJK frame
    num = -sin(phiEclip*d2r)*sin(oblEclip*d2r)+cos(phiEclip*d2r)*cos(oblEclip*d2r)*sin(lonEclip*d2r);
    den = cos(phiEclip*d2r)*cos(lonEclip*d2r);
    AmoonR = atan2(num,den);
    AmoonD = AmoonR/d2r;
else
    DmoonD = [];
    AmoonD = [];
end

end