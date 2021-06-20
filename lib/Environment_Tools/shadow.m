function eclipse = shadow(r,r_sun)
% This function determines if a satellite is in shadow given its position
% vector and the vector to the sun.
% 
% Inputs:
%   r - Vector to spacecraft
%   r_sun - Vector to sun
%   ***Note, vectors can be in any frame as long as they are in the same
%      frame
% 
% Output:
%   illum - returns 0 if satellite illuminated, 1 if satellite in umbral
%   eclipse, 2 if satellite in penumbral eclipse
% 
% Based on Vallado (4th ed) Algorithm 34
% 
% Andrew Dianetti
% 4 March 2019

alpha_umb=0.264121687; %degrees
alpha_pen=0.269007205; %degrees
Re=6378.137; %radius of earth, km

eclipse=0;

u_sat=r/norm(r); %unit vector to satellite
u_sun=r_sun/norm(r_sun); %unit vector to sun

theta_sun_sat=acos(dot(-u_sun,u_sat)); %angle between sun projection vector
% (direction sunlight is shining) and target vector

sat_horiz=norm(r)*cos(theta_sun_sat);
sat_vert=norm(r)*sin(theta_sun_sat);

x=Re/sind(alpha_pen);
pen_vert=tand(alpha_pen)*(x+sat_horiz);

if sat_vert <= pen_vert
    eclipse = 2;
    
    y=Re/sind(alpha_umb);
    umb_vert=tand(alpha_umb)*(y-sat_horiz);
    if sat_vert <= umb_vert
        eclipse = 1;
    end
end