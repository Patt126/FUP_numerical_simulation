function [x,v] = rk4(N, alpha, initial_condition, MS_cost,T,dt)

x = zeros(length(T),N);
x(1,:) = initial_condition(N+1:2*N);
v = zeros(length(T),N);
v(1,:) = initial_condition(1:N);

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
    f1 = F(T(t),x(t,:)');
    f2 = F(T(t)+dt/2,x(t,:)'+dt*f1/2);
    f3 = F(T(t)+dt/2,x(t,:)'+dt*f2/2);
    f4 = F(T(t+1),x(t,:)'+dt*f3);
    v(t+1,:) = v(t,:) + dt*(f1+2*f2+2*f3+f4)'/6;
    g1 = I(T(t),v(t+1,:)');
    g2 = I(T(t)+dt/2,v(t+1,:)'+dt*g1/2);
    g3 = I(T(t)+dt/2,v(t+1,:)'+dt*g2/2);
    g4 = I(T(t+1),v(t+1,:)'+dt*g3);
    x(t+1,:) = x(t,:) + dt*(g1+2*g2+2*g3+g4)'/6;
end
end





