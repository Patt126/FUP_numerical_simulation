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
    T = [0:current_dt:T_max];
    [x,v] = Leapfrog(N, alpha, initial_condition, MS_cost,T,current_dt);
    
        % Compute mode_k, speed_k, energy_k, and total_energy
        mode_k = zeros(length(T), N);  % Initialize mode_k
        speed_k = zeros(length(T), N);  % Initialize speed_k
        energy_k = zeros(length(T), N);  % Initialize energy_k
        total_energy = zeros(length(T), 1);  % Initialize total_energy

    for t = 1:length(T)
        mode_k(t,:) = (A * x(t,:)')';  % Adjust based on your solution matrix
        speed_k(t,:) = (A * v(t,:)')';  % Adjust based on your solution matrix
        energy_k(t,:) = 0.5 * (speed_k(t,:)).^2 + 0.5 * (omega_k(1,:) .* mode_k(t,:)).^2; 
        total_energy(t) = sum(energy_k(t,:));
    end

    formatted_dt = sprintf('%.6f', current_dt);
    result_file = sprintf('results_%s.csv', formatted_dt);
  

    % Read data from CSV file into a matrix
    exact_matrix = readmatrix(result_file);
    solution_matrix = [x(:,1),x(:,5),x(:,16),mode_k(:,1),mode_k(:,3),mode_k(:,5)];

    solutionErrors(solution_matrix,exact_matrix,total_energy(:),1,...
        current_dt,T_max,T,MS_cost);

    % Halve the number of steps for the next iteration
    current_dt = current_dt/2;
end





