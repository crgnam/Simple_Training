function [dX] = collisionDynamics(~, X, g, J, a_applied, L)

% Recover States:
q = X(1:4);
w = X(5:7);
% r = X(8:10);
v = X(11:13);
a = [0; 0; -g];

% Euler's Rotational Equations:
dw = J\(L - cross(w,(J*w)));

% Attitude Kinematics:
Bq = zeros(4,3);
Bq(1:3,:) = cpm(q(1:3)) + diag(repmat(q(4),1,3));
Bq(4,:) = -q(1:3);
dq = (1/2)*Bq*w;

% Form differential state vector:
a = a + a_applied;
dX = [dq; dw; v; a];

end