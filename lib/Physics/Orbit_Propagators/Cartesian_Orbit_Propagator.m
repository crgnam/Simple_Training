% Cartesian Orbit Propagator
% 
% This function replicates the functionality of the previous
% Cartesian_Orbit_Propagator_Script, which was written by Andrew Dianetti 8
% June 2015 to implement functions by Chris Shelton 11 July 2014.  However,
% this updated script is a complete rewrite.
% 
% Andrew Dianetti
% 4 March 2019
% 
% Form:
% [R, V] = Cartesian_Orbit_Propagator(r0,v0,JD_UT1,Density,method,AMR)
% 
% Inputs:
%   r0 - Initial position (column) vector in ECI frame (km)
%   v0 - Initial velocity (column) vector in ECI frame (km)
%   JD_UT1 - Time expressed in Julian Date in UT1
%       Note: Difference between UTC and UT1 can be found here:
%       https://www.nist.gov/pml/time-and-frequency-division/atomic-standards/leap-second-and-ut1-utc-information
%   Density - Density table for drag equations.  Currently using Table 8-4
%       from Vallado 4th edition but can use any table desired.
%   dt - Time step (sec)
%   method - Propagation model:
%       1 - No perturbations (pure Keplarian motion)
%       2 - J2 perturbations
%       3 - J2 + Drag
%       4 - J2-J5 + Sun + Moon + Drag
%   AMR - area-to-mass ratio, m^2/kg
% 
% Outputs:
%   R - propagated ECI position vector (column), km
%   V - propagated ECI velocity vector (column), km

function [R, V]=Cartesian_Orbit_Propagator(r0,v0,JD_UT1,Density,dt,method,AMR)

    Re=6378.137; %Radius of earth
       
    r=norm(r0); %Magnitude of position vector
    
    % Density lookup, using Vallado Eq. 8-33
    % Note: This is fairly crude and could be replaced with a more accurate
    % drag model.  However it should be more accurate than what was in the
    % prior script.
    
    h_ellip=r-Re; %Height of satellite (km)
    ind=max(find(h_ellip>Density(:,1)));
    h_0=Density(ind,1); %Base altitude (km)
    rho_0=Density(ind,2); %Nominal Density (kg/m^3)
    H=Density(ind,3); %Scale Height (km)
%     keyboard
    try
    rho=rho_0*exp(-(h_ellip-h_0)/H); %Density, kg/m^3

    catch
         disp('hi')   
    end

    % Make sure r0 and v0 are column arrays:
    r0=[r0(1); r0(2); r0(3)];
    v0=[v0(1); v0(2); v0(3)];
    
    % RK4 integration:
    dx1=dt*v0;
    dv1=dt*Cartesian(r0,v0,rho,JD_UT1,AMR,method);
    dx2=dt*(v0+dv1/2);
    dv2=dt*Cartesian(r0+dx1/2,v0+dv1/2,rho,JD_UT1,AMR,method);
    dx3=dt*(v0+dv2/2);
    dv3=dt*Cartesian(r0+dx2/2,v0+dv1/2,rho,JD_UT1,AMR,method);
    dx4=dt*(v0+dv3);
    dv4=dt*Cartesian(r0+dx3,v0+dv1,rho,JD_UT1,AMR,method);
    dx=(dx1+2*dx2+2*dx3+dx4)/6;
    dv=(dv1+2*dv2+2*dv3+dv4)/6;

    R=r0+dx;
    V=v0+dv;
    
end

function a = Cartesian(r0,v0,rho,JD_UT1,AMR,method)
    % Define constants for orbit propagation
    mu=398600;
    Re=6378.137; %Radius of earth

    r=norm(r0); %Magnitude of position vector
    
    % Acceleration due to gravity (pure Keplarian motion)
    T1 = -mu/r^3;
    a1 = T1*r0;
    
    if method == 1 %Keplerian only
        a = a1;
    
    else
        % J2 constant
        J2=1.08263e-3;

        % Components of position vector
        x=r0(1);
        y=r0(2);
        z=r0(3);
        
        % J2 acceleration:
        T2=mu*J2*Re^2/2;
        aJ2=[0;0;0]; %Preallocate
        aJ2(1)=T2*(15*z^2/r^7 - 3/r^5)*x;
        aJ2(2)=T2*(15*z^2/r^7 - 3/r^5)*y;
        aJ2(3)=T2*(15*z^2/r^7 - 9/r^5)*z;
        if method == 2 %J2
            a = a1 + aJ2;
        
        else
            % Drag acceleration
            Cd = 2.2; %Drag coefficient.  Per Vallado p. 549, approx 2.2 for
                % satellites in upper atmosphere.  Spheres are around 2.0 to
                % 2.1.
            BC = 1/(AMR*Cd); %Ballistic coefficient (kg/m^2)
            w = 7.292115146706979e-5; %Rotation rate of earth (rad/s), 
                % neglecting Length of Day, which has small impact (on
                % order of 10^-11).  Ref. Vallado pg. 227.
            v_rel = [v0(1)+w*r0(2); v0(2)-w*r0(1); v0(3)]; %Velocity relative 
                % to Earth's atmosphere.  Ref. Vallado p. 550.  Fairly
                % crude approximation, w (rotation rate of atmosphere) will
                % change based on altitude and will be slightly slower than
                % Earth.
            v_rel_m=v_rel*1000; %Relative velocity in m/s (for next equation)
            ad_m = -0.5*1/BC*rho*(v_rel_m'*v_rel_m)*v_rel_m/norm(v_rel_m); %Acceleration due to drag, m/s^2
            ad = ad_m/1000; %Acceleration due to drag in km/s^2
            
            if method == 3 %J2 + Drag (Intended for LEO)
                a = a1 + aJ2 + ad;
            else %method == 4, J2-J5, Sun & Moon gravity, SRP, Drag 
                % (Intended for GEO, but if include drag also good for LEO)
                
                % Constants
                J3=-2.5e-6;
                J4=-1.6e-6;
                J5=-0.15e-6;
                muSun=132712000000;
                muMoon=4903;
    
                % Acceleration due to J3-J5
                % Ref: Bate, Mueller, White
                zR=z/r;
                
                %J3
                aJ3=[0;0;0];
                aJ3(1)=T1*x*J3*2.5*((Re/r)^3)*(3*zR-7*zR^3);
                aJ3(2)=(y/x)*aJ3(1);
                aJ3(3)=T1*z*J3*1.5*((Re/r)^3)*(10*zR-(35/3)*zR)-T1*J3*1.5*((Re/r)^3)*r;
                
                %J4
                aJ4=[0;0;0];
                aJ4(1)=-T1*x*J4*(5/8)*((Re/r)^4)*(3-42*zR^2+63*zR^4);
                aJ4(2)=(y/x)*aJ4(1);
                aJ4(3)=-T1*z*J4*(5/8)*((Re/r)^4)*(15-70*zR^2+63*zR^4);
                
                %J5
                aJ5=[0;0;0];
                aJ5(1)= -T1*x*J5*(3/8)*((Re/r)^5)*(35*zR-210*zR^3+231*zR^5);
                aJ5(2) = (y/x)*aJ5(1);
                aJ5(3) = -T1*z*J5*(1/8)*((Re/r)^5)*(315*zR-945*zR^3+693*zR^5)+T1*J5*(1/8)*((Re/r)^5)*15*r;
                
                % Sun gravity
                r_sun=sunvec(JD_UT1); %Vector to sun
                a_sun1=(r0-r_sun)/norm(r0-r_sun)^3;
                a_sun2=r_sun/norm(r_sun)^3;
                a_sun=-muSun*(a_sun1+a_sun2);
                
                % Moon gravity
                r_moon=moonvec(JD_UT1); %Vector to moon (approximate, assuming TDB=UT1)
                a_moon1=(r0-r_moon)/norm(r0-r_moon)^3;
                a_moon2=r_moon/norm(r_moon)^3;
                a_moon=-muMoon*(a_moon1+a_moon2);
                
                % Solar Radiation Pressure
                eclipse = shadow(r0,r_sun); %returns 0 if illuminated, 1 if in umbra, 2 if in penumbra
                if eclipse ~= 0 %if illuminated
                    Psr = 4.56e-6; %N/m^2
                    cr = 1; %coefficient of reflectivity, this can be changed based on surface model
                    a_SRP=Psr*cr*AMR*(r0-r_sun)/norm(r0-r_sun);
                else
                    a_SRP=[0;0;0];
                end
                
                % Sum all accelerations
                a = a1 + aJ2 + aJ3 + aJ4 + aJ5 + ad + a_sun + a_moon + a_SRP;
                
            end
        end
    end
end