function [x,v] = Nystrom(N, alpha, initial_condition, MS_cost,T,dt)

x = zeros(length(T),N);
x(1,:) = initial_condition(N+1:2*N);
v = zeros(length(T),N);
v(1,:) = initial_condition(1:N);
gamma = (1/12)*(2-nthroot(4,3)-nthroot(16,3));
c = [0,0.5-gamma,0.5,0.5+gamma];
b = [1/(24*gamma^2),1-1/(12*gamma^2),1/(24*gamma^2)];

e = ones(N,1);
linear_part = spdiags([e -2*e e],-1:1,N,N);
next_quadratic = spdiags([e -e], 0:1 ,N,N);
prev_quadratic = spdiags([-e e], -1:0 ,N,N);

F = @(t,state) MS_cost*linear_part*state + alpha*((next_quadratic*state).^2-(prev_quadratic*state).^2);
I = @(t,state) diag(ones(N,1))*state;

f = zeros(4,N);
g = zeros(4,N);
%Apply nystrom scheme
for t = 1:length(T)-1
    f(1,:) = x(t,:);
    g(1,:) = v(t,:);

    for i = 1:3
        f(i+1,:) = f(i,:) + dt*(c(i+1)-c(i))*g(i,:);
        g(i+1,:) = g(i,:) + dt*b(i)*F(t,f(i+1,:)')';
    end
    x(t+1,:) = f(4,:) + dt*(1-c(4))*g(4,:);
    v(t+1,:) = g(4,:);
end



