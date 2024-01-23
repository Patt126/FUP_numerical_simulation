function [x,v] = Ruth3(N, alpha, initial_condition, MS_cost,T,dt)

x = zeros(length(T),N);
x(1,:) = initial_condition(N+1:2*N);
v = zeros(length(T),N);
v(1,:) = initial_condition(1:N);

b_x = [2/3,-2/3,1];
b_v = [7/24,3/4,-1/24];


e = ones(N,1);
linear_part = spdiags([e -2*e e],-1:1,N,N);
next_quadratic = spdiags([e -e], 0:1 ,N,N);
prev_quadratic = spdiags([-e e], -1:0 ,N,N);

F = @(t,state) MS_cost*linear_part*state + alpha*((next_quadratic*state).^2-(prev_quadratic*state).^2);
I = @(t,state) diag(ones(N,1))*state;

f = zeros(4,N);
g = zeros(4,N);

for t = 1:length(T)-1
    g(1,:) = v(t,:);
    f(2,:) = x(t,:); 
   for i =1:3
        g(i+1,:) = g(i,:) +dt*b_v(i)*F(t,f(i+1,:)')';
        f(i+2,:) = f(i+1,:) + dt*b_x(i)*I(t,g(i+1,:)')';
   end
   v(t+1,:) = g(4,:)9hoho;
   x(t+1,:) = f(5,:);

end

end

