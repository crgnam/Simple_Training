function [magField_B,magField_I] = magnetometer(q,sigma,lat,lng,alt,t)
%This function creates measurements for a magnetomer on sat at position pos (lat,lng,alt). Magnetometer readings are created by first generating an magnetic
%field vector in the NED frame using the IGRF model. The magnetic field vector is rotated from the NED frame to the 
%reference frame. The magnetic field vector in the reference frame is stored in magField_I. The magnetic field is rotated
%into the body frame, and after the addition of 3-sigma noise is stored in magField_B.

%Initialize matrices to store inertial and measured magnetic field vectors
magField_I = zeros(3,1);
magField_B = zeros(3,1);

%Get magnetic field vector in NED frame
magField_NED = igrf(t,lat,lng,alt,'geodetic');

%Rotate magnetic field vector from NED to reference (Inertial?) frame
%REPLACE eye(3) w/ actual rotation matrix once reference frame decided
magField_i = eye(3)*magField_NED';

%Generate rotation matrix to take vector from reference to body frame
A_bi = quaternionToAttitudeMatrix(q);

%The magnetic field vector in the body frame (prior to noise)
magField_b = A_bi*magField_i;                                               %ASK CHRIS IF I SHOULD NORM PRIOR TO NOISE ADD!!

%Add noise to magnetic field vector in body frame
z1=(magField_b(1)/abs(magField_b(3)));z2=(magField_b(2)/abs(magField_b(3)));
rz=sigma^2/(1+z1^2+z2^2)*[(1+z1^2)^2 (z1*z2)^2;(z1*z2)^2 (1+z2^2)^2];
[u_z,e_z]=eig(rz);
noise1=sqrt(e_z(1,1))*randn(1);noise2=sqrt(e_z(2,2))*randn(1);
noise_corr=u_z*[noise1;noise2];
alpm=z1+noise_corr(1);betam=z2+noise_corr(2);

%Normalized measured magnetic field vector in body frame
magField_b =[alpm;betam;1*sign(magField_b(3))]/norm([alpm;betam;1]);

%Normalize magneitc field vector in inertial frame
magField_i = magField_i/norm(magField_i);

%Store inertial and measured magnetic field vectors
magField_B(:,1) = magField_b';
magField_I(:,1) = magField_i';

end

