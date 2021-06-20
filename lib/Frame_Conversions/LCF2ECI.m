function [A, ECI] = LCF2ECI(ECI)

    % Convert Cartesian ECI to Spherical:
    [az, el, ~] = cart2sph(ECI(1), ECI(2), ECI(3));
    dec = pi/2 - el;
    
    % Create Upwards Vector:
    U = [cos(az)*sin(dec);
         sin(az)*sin(dec);
         cos(dec)];

    % Create Westward Vector:
    W = -[-sin(az);
           cos(az); 
              0];

    % Create North Vector:
    N = -[cos(az)*cos(dec);
          sin(az)*cos(dec);
             -sin(dec)];

    % Create Local Coordiante Frame:
    A = [N, W, U];

end