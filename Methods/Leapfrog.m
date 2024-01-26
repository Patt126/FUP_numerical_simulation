function [x,v] = Leapfrog(N, alpha, initial_condition, MS_cost,T,dt)



x = zeros(length(T),N);
v = zeros(length(T),N);
v(1,:) = initial_condition(1:N);
x(1,:) = initial_condition(N+1:2*N);
x(2,:) = initial_condition(N+1:2*N) + dt * initial_condition(1:N);


e = ones(N,1);
linear_part = spdiags([e -2*e e],-1:1,N,N); %define equation matrix
next_quadratic = spdiags([e -e], 0:1 ,N,N);
prev_quadratic = spdiags([-e e], -1:0 ,N,N);

F = @(t,x) sqrt(MS_cost)*linear_part*x + alpha*((next_quadratic*x).^2-(prev_quadratic*x).^2);
%
% Solving the equation of motion
%
for t = 3:length(T)
    x(t,:) = 2*x(t-1,:) - x(t-2,:) + (dt^2)*F(t-1,x(t-1,:)')';
end

for t = 2:length(T)
    v(t,:) = (x(t,:)-x(t-1,:))/dt;
end
