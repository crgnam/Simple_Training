function [stars_B,stars_I,numStars] = starTracker(q,FOV,MTH,sigma,maxStars)
%This function creates measurements for a star tracker with given field of view (FOV) and magnitude threshold (MTH) 
%over a specified period of time. The star tracker's optical axis is assumed to be aligned with the satellite body's 
%z-axis. The quaternion representation of the satellite's attitude (w/ respect to the inertial reference frame) at a 
%given timestep is given by q(:,i). Star tracker measurements are computed by adding 3-sigma noise to the "truth". 
%The maximum number of stars capable of being tracked at a single instance is maxStars.

%======= SETUP FOR STAR TRACKER SIMULATION =========================================================================== 

%Determine the location of the star tracker's corners in the sensor frame([x y z]') using its FOV. Corner locations 
%are purposely adjusted to exclude sensor edges and ensure that resulting vectors are of unit length. 

FOV_rad = FOV*(pi/180);                         %Convert the field of view to radians

adj = (2*cos(FOV_rad))/(cos(FOV_rad)+1);        %Adjustment to reduce corner locations to "usable" portion of sensor,
                                                %and ensure unit length.

xDist = sin(FOV_rad/2);                         %Magnitude of distance to sensor corner along the x-axis 
yDist = cos(FOV_rad/2)*sqrt(1-adj);             %Magnitude of distance to sensor corner along the y-axis
zDist = cos(FOV_rad/2)*sqrt(adj);               %Magnitude of distance to sensor corner along the z-axis

c1_s=[-xDist; +yDist; zDist];                   %Location of first corner in the sensor frame                   
c2_s=[+xDist; +yDist; zDist];                   %Location of second corner in the sensor frame
c3_s=[+xDist; -yDist; zDist];                   %Location of third corner in the sensor frame
c4_s=[-xDist; -yDist; zDist];                   %Location of fourth corner in the sensor frame

%Determine the location of the star tracker's corners in the body frame([x y z])

A_bs = eye(3);                                  %The attitude matrix to rotate sensors in the sensor frame to the body
                                                %frame (THIS WILL NEED TO BE UPDATED ONCE POSITION OF SENSOR IS DECIDED)
opticalAxis = A_bs(:,3);                        %A vector representation of the sensor's optical axis

c1_b = A_bs*c1_s;                               %Location of first corner in the body frame
c2_b = A_bs*c2_s;                               %Location of second corner in the body frame
c3_b = A_bs*c3_s;                               %Location of third corner in the body frame
c4_b = A_bs*c4_s;                               %Location of fourth corner in the body frame

%Compute unit normal vectors to the sides of the sensor (in the body frame)

norm_12 = cross(c1_b,c2_b)/sin(FOV_rad);        %unit normal vector for side between corners 1 and 2
norm_23 = cross(c2_b,c3_b)/sin(FOV_rad);        %unit normal vector for side between corners 2 and 3
norm_34 = cross(c3_b,c4_b)/sin(FOV_rad);        %unit normal vector for side between corners 3 and 4
norm_41 = cross(c4_b,c1_b)/sin(FOV_rad);        %unit normal vector for side between corners 4 and 1

%Compute the minimum interstar with respect to the optical axis

minInterstarCosine = c1_b'*c3_b;                        %Minimum interstar cosine (also max since sensor is square)
minInterstar = cos(acos(minInterstarCosine)/2);         %Minimum interstar with respect to the optical axis 

%Read the star map and extract observable stars (based on magnitude threshold)

load mappar;                                            %Load the star map
obsStar_Indices=find(MAGN<=MTH);                        %Find indices of observable stars
ObsStar=VS(:,obsStar_Indices);                          %Extract the observable stars

%======= START STAR TRACKER SIMULATION ==============================================================================

stars_I = zeros(1,3*maxStars);                          %Initialize a matrix that stores the star vectors in the inertial frame
stars_B = zeros(1,3*maxStars);                          %Initialize a matrix that stores the measured star vectors in the body frame

A_bi = quaternionToAttitudeMatrix(q);                   %Get matix representation of the satellite's attitude (transforms 
                                                        %vectors from the inertial frame into the body frame).
numStars = 0;                                           %Initially set the number of stars found by the tracker to 0
ObsStar_Int = (opticalAxis'*A_bi)*ObsStar;              %Calculate interstar w/ respect to optical axis for each star observed                                                   
starIndices = find(ObsStar_Int >= minInterstar);        %Record indices of stars w/ interstar greater than minimum                                                   
    
if (length(starIndices) > 0)                            %If stars w/ interstar greater than the minimum exist
        
    for j = starIndices                                 %Iterate through stars in the list

        star_i = ObsStar(:,j);                          %Get the true star vector in the inertial frame
        star_b = A_bi*star_i;                           %Get the true star vector in the body frame (w/o noise)


        %Generate a vector representation of measured star vector in body frame by adding 3-sigma noise
        z1=(star_b(1)/star_b(3));z2=(star_b(2)/star_b(3));
        rz=sigma^2/(1+z1^2+z2^2)*[(1+z1^2)^2 (z1*z2)^2;(z1*z2)^2 (1+z2^2)^2];
        [u_z,e_z]=eig(rz);
        noise1=sqrt(e_z(1,1))*randn(1);noise2=sqrt(e_z(2,2))*randn(1);
        noise_corr=u_z*[noise1;noise2];
        alpm=z1+noise_corr(1);betam=z2+noise_corr(2);

        star_b =[alpm;betam;1]/norm([alpm;betam;1]);    %Update value to measured star vector in the body frame

        %Check if star is on correct side of sensor plane, and that number of stars being tracked at the given
        %step is less than the max.
        if ( ((norm_12'*star_b) < 0) && ((norm_23'*star_b) < 0) && ((norm_34'*star_b) < 0) && ((norm_41'*star_b) < 0) && (numStars < maxStars) )

            ind_s = 1+(3*numStars);                 %The starting index for storage of the star's information
            stars_B(1,ind_s:ind_s+2) = star_b';     %Store the measured star vector

            star_i = star_i*(1/norm(star_i));       %Make sure the star vector in the inertial frame is unit length
            stars_I(1,ind_s:ind_s+2) = star_i';     %Store the corresponding inertial star vector
            numStars = numStars + 1;                %Increment the number of stored stars (in the timestep) by 1 
        end    
    end        
end
end

