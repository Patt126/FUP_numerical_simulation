function [sol_rk4] = rk4_solver(N, alpha, initial_condition, MS_cost,T,dt)

sol_pos = zeros(length(T),N);
sol_pos(1,:) = initial_condition(N+1:2*N);
sol_vel = zeros(length(T),N);
sol_vel(1,:) = initial_condition(1:N);

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
I = @(t,state) diag(ones(N,1))*state;

for t = 1:length(T)-1
    f1 = F(T(t),sol_pos(t,:)');
    f2 = F(T(t)+dt/2,sol_pos(t,:)'+dt*f1/2);
    f3 = F(T(t)+dt/2,sol_pos(t,:)'+dt*f2/2);
    f4 = F(T(t+1),sol_pos(t,:)'+dt*f3);
    sol_vel(t+1,:) = sol_vel(t,:) + dt*(f1+2*f2+2*f3+f4)'/6;
    g1 = I(T(t),sol_vel(t+1,:)');
    g2 = I(T(t)+dt/2,sol_vel(t+1,:)'+dt*g1/2);
    g3 = I(T(t)+dt/2,sol_vel(t+1,:)'+dt*g2/2);
    g4 = I(T(t+1),sol_vel(t+1,:)'+dt*g3);
    sol_pos(t+1,:) = sol_pos(t,:) + dt*(g1+2*g2+2*g3+g4)'/6;
end
sol_rk4 = [sol_vel,sol_pos];
end





