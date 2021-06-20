function[X0] = elementToState(inc,e,v,RAAN,w,a,mu,~)
%% Calculate Eccentric anomaly:
% Use Newton-Raphson method to solve for E
f = @(E) E-e*sin(E)-v;
df = @(E) 1-e*cos(E);

%Initialize values:
E = 1;
F = 1;
c = 1;
eps = .001;

while abs(F) > eps && abs(c) > eps
    F = f(E);
    dF = df(E);
    c = -F/dF;
    E = E + c;
end

%% Prep the positions and velocities:
rmag = a*(1-e*cos(E));

x = a*(cos(E)-e);
y = a*sqrt(1-e^2)*sin(E);
xdot = -sqrt(mu*a)/rmag*sin(E);
vdot = sqrt(mu*a*(1-e^2))/rmag*cos(E);

%% Transformation from equatorial coordinate system to perifocal:
sRAAN = sin(RAAN);
cRAAN = cos(RAAN);

sw = sin(w);
cw = cos(w);
sinc = sin(inc);
cinc = cos(inc);

Q11 = cRAAN*cw - sRAAN*cinc*sw; 
Q12 = sRAAN*cw + cRAAN*sw*cinc; 
Q13 = sw*sinc;

Q21 = -cRAAN*sw - sRAAN*cw*cinc;
Q22 = -sRAAN*sw + cRAAN*cw*cinc;
Q23 = cw*sinc;

Q31 = sRAAN*sinc;
Q32 = -cRAAN*sinc;
Q33 = cinc;

Q = [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];

%% Produce the initial state vectors:
r0 = Q'*[x y 0]';
v0 = Q'*[xdot vdot 0]';
X0 = [r0; v0];
    
end