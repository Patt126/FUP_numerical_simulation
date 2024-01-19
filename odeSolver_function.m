function [T_sol,sol_ode15s] = odeSolver(N,alpha,initial_condition ,MS_cost,T)

sol_ode15s = zeros(length(T),2*N); % In a row first N  are position, last N are speed

%Define  matrix for the system
linear_part = zeros(2*N,2*N);
next_quadratic = zeros(2*N,2*N);
prev_quadratic = zeros(2*N,2*N);

for i = 1:N
    linear_part(i,i+N) = -2*MS_cost;
    next_quadratic(i,i+N) = 1;
    prev_quadratic(i,i+N) = 1;
    if i ~= 1
        linear_part(i,i-1+N) = 1*MS_cost;
        prev_quadratic(i,i+N-1) = -1;
    end
    if i ~= N
        linear_part(i,i+1+N) = 1*MS_cost;
        next_quadratic(i,i+N+1) = -1;
    end
    linear_part(i+N,i) = 1;
end

option = odeset('RelTol',1e-4,'AbsTol',1e-6);
F = @(t,state) linear_part*state + alpha*((next_quadratic*state).^2-(prev_quadratic*state).^2);
[T_sol,sol_ode15s] = ode15s(F,T,initial_condition,option);

end

