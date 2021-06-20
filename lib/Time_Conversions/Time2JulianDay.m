function [ JD ] = Time2JulianDay( time )
% This function converts a six element time array consisting of:
%   year, month, day, hour, minute, second
% into the Julian Date, JD. This date is used for many calculations in
% astronomy.

% First unpack the time array into its components for the sake of clarity
year = time(1);
month = time(2);
day = time(3);
hour = time(4);
minute = time(5);
sec = time(6);

% Now break the conversion up into different terms 
term2 = floor( 7*(year + floor( (month+9) /12) ) /4 );
term3 = floor( 275*month / 9);

% Perform the calculation
JD = 367*year - term2 + term3 + day + 1721013.5 + hour/24 + minute/1440 + sec/86400;

end