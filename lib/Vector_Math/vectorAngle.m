% This function returns the angle, in degrees, that exists between two
% vectors that are expressed in the same reference frame.
% Written by Nick DiGregorio

function theta = vectorAngle(v1, v2)

theta = acosd( dot(v1,v2) / ( norm(v1)*norm(v2) ) );

end