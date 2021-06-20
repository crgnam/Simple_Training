% This function is not used in the Master Sim, but it IS part of the mission
% module, because it's easier to convert this than to rewrite it in C++.
function [sunPos,sunVel,sunVec,earthVec,sunVecBody] = CalcSunEarthVectors(JD,pos,q)
    [sunPos, sunVel] = CalculateSunRV(JD);
    sunPos = sunPos';
    sunVel = sunVel';
    
    sunVec = sunPos - pos;
    %sunVec = q2a(q)*sunVec;
    sunVec = sunVec/norm(sunVec);
    
    earthVec = -pos;
    %earthVec = q2a(q)*earthVec;
    earthVec = earthVec/norm(earthVec);
    
    sunVecBody = q2a(q)*sunVec;
end
