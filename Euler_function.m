function [sol_Euler] = EulerMethod(N, alpha, initial_condition, MS_cost,T,dt)

sol_Euler = zeros(length(T),2*N);
sol_Euler(1,:) = initial_condition;

linear_part = zeros(N,N);
next_quadratic = zeros(N,N);
prev_quadratic = zeros(N,N);
for i = 1:N
    linear_part(i,i) = -2;
    next_quadratic(i,i) = 1;
    prev_quadratic(i,i) = 1;
    if i ~= 1
        linear_part(i,i-1) = 1;
        prev_quadratic(i,i-1) = -1;
    end
    if i ~= N
        linear_part(i,i+1) = 1;
        next_quadratic(i,i+1) = -1;
    end
end

F = @(t,state) linear_part*state + alpha*((next_quadratic*state).^2-(prev_quadratic*state).^2);

for t = 2:length(T)
    sol_Euler(t,1:N) = sol_Euler(t-1,1:N) + dt*F(t,sol_Euler(t-1,N+1:2*N)')';
    sol_Euler(t,N+1:2*N) = sol_Euler(t-1,N+1:2*N) + dt*sol_Euler(t,1:N);
end

