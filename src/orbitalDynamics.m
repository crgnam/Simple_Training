function [dX] = orbitalDynamics(~,X,mu)
    % States:    
    r = X(1:3);
    v = X(4:6);
    
    % Accelerations:
%     a_drag = dragModel(r,v);
    a = -(mu/norm(r)^3)*r;
%     a_j2 = j2Model(r);
%     a = a_drag + a_grav + a_j2;
    
    % Differential states:
    dX = [v; a];
end