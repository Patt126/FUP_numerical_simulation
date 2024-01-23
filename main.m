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
omega_0 = 2*sin(pi/(2*(N+1))); %consider the possibility of compute here normal mode and then pass to function
A = zeros(N,N);
for row = 1:N
    A(row,:) = sqrt(2/(N+1))*sin((pi*row.*(1:N))/(N+1)); 
end
mode_0(1) = 1*sqrt(2)/omega_0;
initial_condition(N+1:2*N,1) = A\mode_0;
%

% [T_sol,sol_ode15s] = odeSolver(N, alpha, initial_condition, MS_cost, T); 
% [energy_k,total_energy_error] = plotEnergy(sol_ode15s(:,N+1:2*N),sol_ode15s(:,1:N),T_sol,N,A,n_mode,'Ode15s');
% disp("Ode15: ")
% disp(max(total_energy_error))
% 

%Euler Method
%
[x_Euler,v_Euler] = EulerMethod(N, alpha, initial_condition, MS_cost,T,dt);
[energy_k,total_energy_error] = plotEnergy(x_Euler,v_Euler,T,N,A,n_mode,'Euler');
disp("Euler: ")
disp(max(total_energy_error))

% 
% % Runge-Kutta 3
% [x_rk3,v_rk3] = rk3(N, alpha, initial_condition, MS_cost,T,dt);
% [energy_k,total_energy_error] = plotEnergy(x_rk3,v_rk3,T,N,A,n_mode,'RK-3');
% disp("RK-3: ")
% disp(max(total_energy_error))
% 
% % Runge Kutta 4
% %
% [x_rk4,v_rk4] = rk4(N, alpha, initial_condition, MS_cost,T,dt);
% [energy_k,total_energy_error] = plotEnergy(x_rk4,v_rk4,T,N,A,n_mode,'RK-4');
% disp("RK-4: ")
% disp(max(total_energy_error))

% [x_ruth3,v_ruth3,T_sol] = Ruth3(N, alpha, initial_condition, MS_cost,T,1e-6);
% [energy_k,total_energy_error] = plotEnergy(x_ruth3,v_ruth3,T_sol,N,A,n_mode,'Ruth-3');
% disp("Ruth-3: ")
% disp(max(total_energy_error))



