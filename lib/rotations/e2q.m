function [q, A] = e2q(roll,pitch,yaw)
% This function calculates the shuster quaternion attitude of the input
% roll, pitch yaw euler angles

    if roll == 180
        roll = 179.99;
    end
    if pitch == 180
        pitch = 179.99;
    end
    if yaw == 180
        yaw = 179.99;
    end

    T1 = [cosd(roll)  -sind(roll)  0;
          sind(roll)   cosd(roll)  0;
               0           0       1];
          
    T2 = [1       0           0;
          0  cosd(pitch)  sind(pitch);
          0 -sind(pitch)  cosd(pitch)];

    T3 = [ cosd(yaw) 0  sind(yaw);
               0     1      0;
          -sind(yaw) 0  cosd(yaw)];

    A = T3*T2*T1;
    q = a2q(A);
end