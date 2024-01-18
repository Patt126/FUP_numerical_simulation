function [energy_k,x,T_sol] = odeSolver(N,alpha,n_mode,MS_cost, t_max,number_step)

sqt_MS = sqrt(MS_cost);
dt = t_max/number_step;
T = (0:dt:t_max);

solution = zeros(length(T),2*N); % In a row first N  are position, last N are speed
mode_0 = zeros(N,1); %Initial displacement in normal modes representation

% define normal modes frequencies, coordinates and energy per mode and in
%time
omega_k = 2*sqt_MS*sin(pi.*(1:n_mode)./(2*(N+1))); %frequencies
mode_k = zeros(length(T),n_mode); % modes
energy_k = zeros(length(T),n_mode); %modes' energy

%
% Define the matrix to pass to from cartesian coordinates to normal modes

A = zeros(N,N);
for row = 1:N
    A(row,:) = sqrt(2/(N+1))*sin((pi*row.*(1:N))/(N+1)); 
end
mode_0(1) = 1*sqrt(2)/omega_k(1);
x_0 = A\mode_0;
initial_condition = zeros(2*N,1);
initial_condition(N+1:2*N,1) = x_0;

for i = 1:n_mode
    mode_k(1,i) = sqrt(2/(N+1))*sum(x_0'.*sin(pi*i*(1:N)/(N+1)));
    energy_k(1,i) = 0.5*(omega_k(i)*mode_k(1,i))^2;
end


%Define  matrix for the system
linear_part = zeros(2*N,2*N);
next_quadratic = zeros(2*N,2*N);
prev_quadratic = zeros(2*N,2*N);
for i = 1:N
    linear_part(i,i+N) = -2*sqt_MS;
    next_quadratic(i,i+N) = 1;
    prev_quadratic(i,i+N) = 1;
    if i ~= 1
        linear_part(i,i-1+N) = 1*sqt_MS;
        prev_quadratic(i,i+N-1) = -1;
    end
    if i ~= N
        linear_part(i,i+1+N) = 1*sqt_MS;
        next_quadratic(i,i+N+1) = -1;
    end
    linear_part(i+N,i) = 1;
end

option = odeset('RelTol',1e-4,'AbsTol',1e-6);
F = @(t,state) linear_part*state + alpha*((next_quadratic*state).^2-(prev_quadratic*state).^2);
[T_sol,solution] = ode15s(F,T,initial_condition,option);

total_energy = zeros(length(T),1);
total_energy(1) = sum(energy_k(1,:));
for t = 2:length(T)-1
    for k = 1:n_mode
    mode_k(t,k) = sqrt(2/(N+1))*sum(.*sin(pi*k.*(1:N)/(N+1)));
    energy_k(t,k) = 0.5*((mode_k(t,k)-mode_k(t-1,k))/dt)^2 + 0.5*(omega_k(k)*mode_k(t,k))^2;
    end
    total_energy(t) = sum(energy_k(t,:));
end

x = solution(:,N+1:2*N);

