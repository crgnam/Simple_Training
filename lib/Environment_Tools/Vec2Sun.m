function [RsunKM,DsunD,AsunD] = Vec2Sun(TUT1,ADon)
%   using the TLE year and day calculate the right ascension, declination, and
%   distance vector to the sun from the center of the Earth's Geocentric
%   Equatorial System(IJK) pg. 159 in Vallado
%   this script follows Vallado's Algorithm 29 pg 281
%
%Inputs:
%   TUT1 : present time
%   ADon : boolean (=1 if want sun angles out, else 0)

d2r = pi/180;
AU = 149597870;%1 AU = .. km

%mean longitude of sun in degrees then mod between 0 and 360
mLsD = 280.46 + 36000.771*TUT1;
mLsD = mod(mLsD,360);

%mean anomaly of the sun
mAsD = 357.5277233+35999.05034*TUT1;
mAsD = mod(mAsD,360);

%earth's eccentricity calc
eE = 0.016708617-0.000042037*TUT1-0.0000001236*TUT1^2;

lonEclip = mLsD + 1.914666471*sin(mAsD*d2r) + 0.019994643*sin(2*mAsD*d2r);

oblEclip = 23.439291 - 0.0130042*TUT1;

rSun = 1.000140612 - 0.016708617*cos(mAsD*d2r) - 0.000139589*cos(2*mAsD*d2r);

%distance vector to Sun in AU
%
RsunAU =[rSun*cos(lonEclip*d2r);rSun*cos(oblEclip*d2r)*sin(lonEclip*d2r);rSun*sin(oblEclip*d2r)*sin(lonEclip*d2r)];

%convert astronomical units to kilometers
RsunKM = RsunAU.*AU;

if ADon
    %calculate declination of sun in IJK frame
    phiEclip = 0;
    DsunR = asin(sin(phiEclip*d2r)*cos(oblEclip*d2r)+cos(phiEclip*d2r)*sin(oblEclip*d2r)*sin(lonEclip*d2r));
    DsunD = DsunR/d2r;

    %calculate right ascension of sun in IJK frame
    num = -sin(phiEclip*d2r)*sin(oblEclip*d2r)+cos(phiEclip*d2r)*cos(oblEclip*d2r)*sin(lonEclip*d2r);
    den = cos(phiEclip)*cos(lonEclip*d2r);
    AsunR = atan2(num,den);
    AsunD = AsunR/d2r;
else
    DsunD = [];
    AsunD = [];
end

end
