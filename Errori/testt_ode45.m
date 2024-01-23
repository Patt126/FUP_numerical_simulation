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

lattice_parameter = 1; %lattice parameter
initial_condition = zeros(2*N,1);
%
omega_k = 2*sin(pi.*(1:N)./(2*(N+1))); %frequencies

mode_0 = zeros(N,1);
omega_0 = 2*sin(pi/(2*(N+1))); %consider the possibility of compute here normal mode and then pass to function
A = zeros(N,N);
for row = 1:N
    A(row,:) = sqrt(2/(N+1))*sin((pi*row.*(1:N))/(N+1)); 
end
mode_0(1) = 1*sqrt(2)/omega_0;
initial_condition(N+1:2*N,1) = A\mode_0;
%

T_max = 1e4;
dt_max = 0.5;
dt_min = 7e-3;  % Adjust this value as needed

current_dt = dt_max;

while current_dt >= dt_min
    % Perform your computations with current_dt
    [T_sol, sol_ode45] = odeSolver(N, alpha, initial_condition, MS_cost, 0:current_dt:T_max);

    % Compute mode_k, speed_k, energy_k, and total_energy
    mode_k = zeros(length(T_sol), N);  % Initialize mode_k
    speed_k = zeros(length(T_sol), N);  % Initialize speed_k
    energy_k = zeros(length(T_sol), N);  % Initialize energy_k
    total_energy = zeros(length(T_sol), 1);  % Initialize total_energy

    for t = 1:length(T_sol)
        mode_k(t,:) = (A * sol_ode45(t, N+1:2*N)')';  % Adjust based on your solution matrix
        speed_k(t,:) = (A * sol_ode45(t, 1:N)')';  % Adjust based on your solution matrix
        energy_k(t,:) = 0.5 * (speed_k(t,:)).^2 + 0.5 * (omega_k(1,:) .* mode_k(t,:)).^2; 
        total_energy(t) = sum(energy_k(t,:));
    end

    % Write data to CSV file
    result_data = [sol_ode45(:,33),sol_ode45(:,37),sol_ode45(:,48), mode_k(:,1),mode_k(:,3),mode_k(:,5), total_energy];
    result_file = sprintf('results_%.6f.csv', current_dt);  % Create file name with dt
    writematrix(result_data, result_file);

    % Halve the number of steps for the next iteration
    current_dt = current_dt/2;
end





