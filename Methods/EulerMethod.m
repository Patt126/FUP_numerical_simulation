function [x,v] = EulerMethod(N, alpha, initial_condition, MS_cost, T, dt)

x = zeros(length(T),N);
v = zeros(length(T),N);
v(1,:) = initial_condition(1:N);
x(1,:) = initial_condition(N+1:2*N);

e = ones(N,1);
linear_part = spdiags([e -2*e e], -1:1, N, N);
next_quadratic = spdiags([e -e], 0:1 , N, N);
prev_quadratic = spdiags([-e e], -1:0 , N, N);


F = @(t,x) MS_cost * linear_part * x + alpha * ((next_quadratic*x).^2-(prev_quadratic*x).^2);

for t = 2:length(T)
    v(t,:) = v(t-1,:) + dt * F(t-1,x(t-1,:)')';
    x(t,:) = x(t-1,:) + dt * v(t,:);
end

