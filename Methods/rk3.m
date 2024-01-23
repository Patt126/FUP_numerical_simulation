function [x,v] = rk3(N, alpha, initial_condition, MS_cost,T,dt)

x = zeros(length(T),N);
x(1,:) = initial_condition(N+1:2*N);
v = zeros(length(T),N);
v(1,:) = initial_condition(1:N);

e = ones(N,1);
linear_part = spdiags([e -2*e e],-1:1,N,N);
next_quadratic = spdiags([e -e], 0:1 ,N,N);
prev_quadratic = spdiags([-e e], -1:0 ,N,N);

F = @(t,state) sqrt(MS_cost)*linear_part*state + alpha*((next_quadratic*state).^2-(prev_quadratic*state).^2);
I = @(t,state) diag(ones(N,1))*state;

for t = 1:length(T)-1
   f1 = F(T(t),x(t,:)');
   f2 = F(T(t)+dt,x(t,:)' + dt*f1);
   f3 = F(T(t)+0.5*dt,x(t,:)' +0.25*dt*(f1+f2));
   v(t+1,:) = v(t,:) + dt*(f1+f2+4*f3)'/6;
   g1 = I(T(t),v(t+1,:)');
   g2 = I(T(t)+dt,v(t+1,:)'+dt*g1);
   g3 = I(T(t)+0.5*dt,v(t+1,:)'+0.25*dt*(g1+g2));
   x(t+1,:) = x(t,:) + dt*(g1+g2+4*g3)'/6;
end

end

