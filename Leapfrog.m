clear all
close all
format long
%
% FPUT simulation
%
disp("WELCOME TO FPUT SIMULATION")
disp("For the following total energy, mass and spring constant are suppose to be unitary")
N = input("Number of particles in motion: ");
flag = true;
while flag
alpha = input("Quadratic nonlinearity coefficent: ");
if alpha > 0.5
    disp("Aplha must be a small perturbation. Try with a smaller number.");
else flag = false;
end
end
n_mode = input("Number of modes analyzed: ");
while n_mode > N || n_mode <= 0 
    disp("Number of modes must be a positive integer less or equal to the number of particles in motion");
    n_mode = input("Number of modes analyzed: ");
end
t_max = input("Simulation time: ");
dt = input("time step: ");

lattice_parameter = 1; %lattice parameter
T = (0:dt:t_max);
x = zeros(length(T),N+2); % Displacements from equilibrium position
mode_0 = zeros(N,1); %Initial displacement in normal modes representation
v = zeros(length(T),N+2); % Speed
a = zeros(length(T),N+2); % Acceleration
%
% define normal modes frequencies, coordinates and energy per mode and in
%time
omega_k = 2*sin(pi.*(1:n_mode)./(2*(N+1))); %frequencies
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
x(1,2:N+1) = x_0';
x(2,2:N+1) = x(1,2:N+1);
for i = 1:n_mode
    mode_k(1,i) = sqrt(2/(N+1))*sum(x(1,2:N+1).*sin(pi*i*(1:N)/(N+1)));
    energy_k(1,i) = 0.5*(omega_k(i)*mode_k(1,i))^2;
end
%
% Solving the equation of motion
%
for t = 3:length(T)
    for n = 1:N
      x(t,n+1) = 2*x(t-1,n+1) - x(t-2,n+1) + dt*dt*((x(t-1,n+2)-2*x(t-1,n+1)+x(t-1,n)) + alpha *((x(t-1,n+2)-x(t-1,n+1))^2-(x(t-1,n+1)-x(t-1,n))^2));
    end
end
for t = 2:length(T)
for k = 1:n_mode
    mode_k(t,k) = sqrt(2/(N+1))*sum(x(t,2:N+1).*sin(pi*k.*(1:N)/(N+1)));
    energy_k(t,k) = 0.5*((mode_k(t,k)-mode_k(t-1,k))/dt)^2 + 0.5*(omega_k(k)*mode_k(t,k))^2;
end
end
%
% Calculating the energy at time t
%

%{
for t = 2:length(T)-1
    for k = 1:n_mode
    mode_k(t,k) = sqrt(2/(N+1))*sum(x(t,2:N+1).*sin(pi*k.*(1:N)/N+1));
    energy_k(t,k) = 0.5*(mode_k(t+1,k)-mode_k(t-1,k)/(2*dt))^2 + (omega_k(k)*mode_k(t,k))^2;
    end
end

for i=1:n_mode
    figure(1)
    plot(T,energy_k(:,i),'b-','LineWidth',1);
    hold on;
    xlabel('t');
    ylabel('E_{k}');
    string = strcat('N=',num2str(N),', alpha=',num2str(alpha));
end
%}
plot(T,energy_k(:,1),'r-','LineWidth',1);hold on
plot(T,energy_k(:,2),'k-','LineWidth',1);
plot(T,energy_k(:,3),'b-','LineWidth',1);
plot(T,energy_k(:,4),'m-','LineWidth',1);
plot(T,energy_k(:,5),'g-','LineWidth',1);
xlabel('t');ylabel('E_{k}');
legend('mode 1','mode 2','mode 3','mode 4','mode 5')
%
% Graphic simulation of particles
%
delta = lattice_parameter*ones(length(T),1);
figure(2)
for i = 1:N+1
    x(:,i) = x(:,i) + i*delta;
end
for j=1:length(T)
    figure(2)
    plot(x(j,:),zeros(N+2,1),'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    drawnow;
end






  
   
