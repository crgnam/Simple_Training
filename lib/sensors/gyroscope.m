function [omega_m,bias_m] = gyroscope(omega_t,bias_t,sigma_u,sigma_v,numSteps,dt)

num_g=dt*[1 1];den_g=2*[1 -1];
[phi_g,gam_g,c_g,d_g]=tf2ss(num_g,den_g);
bias1=dlsim(phi_g,gam_g,c_g,d_g,sigma_u/sqrt(dt)*randn(numSteps,1),bias_t*pi/180/3600/dt);
bias2=dlsim(phi_g,gam_g,c_g,d_g,sigma_u/sqrt(dt)*randn(numSteps,1),bias_t*pi/180/3600/dt);
bias3=dlsim(phi_g,gam_g,c_g,d_g,sigma_u/sqrt(dt)*randn(numSteps,1),bias_t*pi/180/3600/dt);

bias_m = [bias1 bias2 bias3];                                                                 %A matrix that stores the true bias at each timestep
omega_m = omega_t+sqrt(sigma_v^2/dt+1/12*sigma_u^2*dt)*randn(numSteps,3)+bias_m;              %A matrix that stores the measured angular velocity at each timestep.

end

