clear all
close all
format long
%
% Main script: FPUT simulation. Different method are presented below:
% Euler method, Leapfrog, Nystrom 3, Runge-Kutta 4, Ruth 4
% The script returns two plot: 
% - the first shows energy trend over time of the first n_modes (you can modify n_modes at line 14) 
% - the second show equipartition
% Moreover is computed the error for a specific coordinate and for the
% total energy, comparing the method with the solution of ode45.
% 
N = 32; % Number of particles
alpha = 0.1; % non-linear coefficient
MS_cost = 1; % Mass-Spring constant
energy = 1; % System's energy
n_mode = 5;  % Number of modes to visualize in the plot
t_max = 1e4; % Simulation time
number_step = 1e6; % Number of step
dt = t_max/number_step; % Step size
T = (0:dt:t_max); % Create time vector
lattice_parameter = 1; % Distance between atoms
initial_condition = zeros(2*N,1); % First N's are velocities, second N's are positions
mode_0 = zeros(N,1); % initialize vector for the first mode
omega_0 = 2*sin(pi/(2*(N+1))); % frequency of the first mode
%
% Create matrix A to pass from cartesian coordinates to nomal modes
%
A = zeros(N,N); 
for row = 1:N
    A(row,:) = sqrt(2/(N+1)) * sin((pi * row .* (1:N))/(N+1)); 
end
mode_0(1) = energy * sqrt(2)/omega_0; % At t=0 energy is all in the first mode --> energy = 0.5*(omega_0*mode_0)^2. Solve for mode_0.
initial_condition(N+1:2*N,1) = A\mode_0; % Assign initial positions
%
%Euler Method
%
[x_Euler, v_Euler] = EulerMethod(N, alpha, initial_condition, MS_cost, T, dt); % Solve with Euler Method
[energy_k, total_energy] = plotEnergy(x_Euler, v_Euler, T, N, A, MS_cost, n_mode, 'Euler'); % Plot Energy modes over time and show equipartition.
[rel_err_x_inf, rel_err_energy_inf] = Error_method(x_Euler(:,14), total_energy, N, alpha, initial_condition , MS_cost, T); % Calculate relative error for a single coordinate and total energy with norm inf.
disp("Relative error with inf-norm for coordinate 14")
disp(rel_err_x_inf);
disp("Relative error with inf-norm for total energy ")
disp(rel_err_energy_inf);

%
% Leapfrog
%
% [x_Leapfrog, v_Leapfrog] = Leapfrog(N, alpha, initial_condition, MS_cost, T, dt); 
% [energy_k, total_energy] = plotEnergy(x_Leapfrog, v_Leapfrog, T, N, A, n_mode, 'Leapfrog'); 
% [rel_err_x_inf, rel_err_energy_inf] = Error_method(x_Leapfrog(:,14), total_energy); 
% 
% Nystrom 3
%
% [x_Nystrom, v_Nystrom] = Nystrom(N, alpha, initial_condition, MS_cost, T, dt); 
% [energy_k, total_energy] = plotEnergy(x_Nystrom, v_Nystrom, T, N, A, n_mode, 'Nystrom'); 
% [rel_err_x_inf, rel_err_energy_inf] = Error_method(x_Nystrom(:,14), total_energy); 
% 
% Runge Kutta 4
% 
% [x_rk4,v_rk4] = rk4(N, alpha, initial_condition, MS_cost, T, dt);
% [energy_k, total_energy] = plotEnergy(x_rk4, v_rk4, T, N, A, n_mode, 'Runge-Kutta 4'); 
% [rel_err_x_inf, rel_err_energy_inf] = Error_method(x_rk4(:,14), total_energy); 
%
% Ruth 3
%
% [x_ruth3,v_ruth3,T_sol] = Ruth3(N, alpha, initial_condition, MS_cost, T, dt);
% [energy_k, total_energy] = plotEnergy(x_ruth3, v_ruth3, T, N, A, n_mode, 'Ruth 3'); 
% [rel_err_x_inf, rel_err_energy_inf] = Error_method(x_ruth3(:,14), total_energy); 
% 




