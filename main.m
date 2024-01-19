clear all
close all
format long
%
% Main script: FPUT simulation
%
N = 32;
alpha = 0.1;
MS_cost = 1; % Mass-Spring constant
n_mode = 5;
t_max = 1e4;
number_step = 1e5;
dt = t_max/number_step;
T = (0:dt:t_max);
lattice_parameter = 1; %lattice parameter
initial_condition = zeros(2*N,1);
%
mode_0 = zeros(N,1);
omega_0 = 2*sin(pi/(2*(N+1)));
A = zeros(N,N);
for row = 1:N
    A(row,:) = sqrt(2/(N+1))*sin((pi*row.*(1:N))/(N+1)); 
end
mode_0(1) = 1*sqrt(2)/omega_0;
initial_condition(N+1:2*N,1) = A\mode_0;
%
% Ode15s Solver
%
[T_sol,sol_ode15s] = odeSolver(N, alpha, initial_condition, MS_cost, T); 
plotEnergy(sol_ode15s,T_sol,N,n_mode);
figure(1)
title('Total Energy trend with ode15s solver');
figure(2)
title('Energy modes over time with ode15s');
%
% Euler Method
%
[sol_Euler] = EulerMethod(N, alpha, initial_condition, MS_cost,T,dt);
plotEnergy(sol_Euler,T,N,n_mode);
figure(3)
title('Total Energy trend with Euler method');
figure(4)
title('Energy modes over time with Euler method');
%
% Runge-Kutta 3
%
[sol_rk3] = rk3_solver(N, alpha, initial_condition, MS_cost,T,dt);
plotEnergy(sol_rk3,T,N,n_mode);
figure(5)
title('Total Energy trend with Runge-Kutta-3');
figure(6)
title('Energy modes over time with Runge-Kutta-3');
%
% Runge Kutta 4
%
[sol_rk4] = rk4_solver(N, alpha, initial_condition, MS_cost,T,dt);
plotEnergy(sol_rk3,T,N,n_mode);
figure(7)
title('Total Energy trend with Runge-Kutta-4');
figure(8)
title('Energy modes over time with Runge-Kutta-4');


