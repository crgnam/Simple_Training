function Bvec = elerot(Nvec,Axis,Angle)
%% elerot = Elementary Rotation
% function will take a vector in a given frame and output the same vector
% after a specified rotation has occured
%
%
% inputs:
% Nvec = 3x1 vector in initial frame
% Axis = number of axis about which rotation is occuring (1,2,3)
% Angle = angle through which rotation occurs in DEGREES
% 
%  outputs:
%  Bvec = 3x1 vector expressed in initial frame coordinates after rotation

% revisions:
%  9/21/15 - initial creation - Nicholas Davidson

Bvec = zeros(3,1);

%% Error messages for improper inputs

% check that the correct amount of input arguements are present
if nargin < 3
    error('not enough input arguements')
end

% check that the axis arguement is compatible (1, 2 or 3)
if Axis ~= 1 && Axis ~= 2 && Axis ~= 3
    error('axis of rotation must be either 1, 2 or 3')
end

% determine size of input Nvec
[m,n] = size(Nvec);

% make sure input vector is a column vector
if n > m
    Nvec = Nvec';
    [m,n] = size(Nvec);
end

% check that Nvec is a single column
if n ~= 1
    error('input must be a vector')
end

% check that Nvec has three components
if m ~= 3
    error('input vector must have three components')
end

%% Carry out specified rotation

% use switch statement to determine which type of rotation was asked for
%
% first initialize the correct transformation matrix and multiply by given
% vector to find needed rotated vector
switch Axis
    case 1
        
        C = [1 0 0; ...
            0 cosd(Angle) -sind(Angle); ...
            0 sind(Angle) cosd(Angle)];
        
        Bvec = C*Nvec;  
        
    case 2
        
        C = [cosd(Angle) 0 sind(Angle); ...
            0 1 0; ...
            -sind(Angle) 0 cosd(Angle)];
        
        Bvec = C*Nvec;
        
    case 3
        
        C = [cosd(Angle) -sind(Angle) 0; ...
            sind(Angle) cosd(Angle) 0; ...
            0 0 1];
        
        Bvec = C*Nvec;
end

end