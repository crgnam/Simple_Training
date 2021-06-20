function E = Eanomallyfunc(e,Me)
%Created 7/11/2014 by Chris Shelton
%% numericaly solve E = Me + e*sin(E)
% This Code uses the newton raphson method to solve for the Eccentric
% anomaly

%Make and Initial Guess
Me_old=Me;
Me = mod(Me,(2*pi));
E = Me;

old = E;
new = 0;
a = old;
b = new;
n = 0;
%Solve for E
while and(abs(b-a) > 1E-13,n < 300) 
    new = old - ((old- Me - (e*sin(old)))/(1-e*cos(old)));
    a = old;
    b = new;
    old=new; 
    n = n+1;
end

%Echo to user if solve doesn't converge
if n == 300
    True_Anomally = 'not convergent'
end
E = new;


