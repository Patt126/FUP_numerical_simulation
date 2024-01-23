
% Example usage:
N = 32;  % Set your value of N
alpha = 0.1;
initial_condition = zeros(N, 1);  % Replace with your initial condition
MS_cost = 1.0;
T = 0:dt:1000;  % Replace with your time array
dt = 0.1;
h0 = 0.1;
hmin = 1e-5;
TOL = 1e-6;

% Call the mirkf45 method
[t, u] = mirkf45(N, alpha, initial_condition, MS_cost, T, dt, h0, hmin, TOL);

