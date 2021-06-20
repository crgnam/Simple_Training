function [directions] = unitDirections(varargin)
    directions = zeros(3,nargin);
    for ii = 1:nargin
        azEl = deg2rad(varargin{ii});
        [x,y,z] = sph2cart(azEl(1), azEl(2), 1);
        directions(:,ii) = [x; y; z];
    end
end