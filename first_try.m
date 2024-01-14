function [Ek, wk, X] = FPUT_simulation(N,alpha,beta,scenario,t_max,dt,plot_opt)
% Inputs: N: number of particles
%         alpha: constant parameter for quadratic nonlinearity model
%         beta: constant parameter for cubic nonlinearity model
%         scenario: 1, %100 energy in mode 1; 2, 75% energy in mode 1 and 25% in mode 3 
%         t_max: maximum simulation time
%         dt: time step, optional, default is 0.1
%         plot_opt: plot optional logical flag, default is 0 (False)
% Outputs: Ek: Energy for each normal modes
%          wk: normal mode freq
%          X: particle positions
% example: [Ek,wk,X] = FPUT_simulation(32,0,3,2,8000,0.1,1);

clc;
a = 10;
% define time parameters
if nargin <= 5
   dt = 0.1; % time step, defaul time step 0.1s 
end
T = (0:dt:t_max); % time span

if nargin <= 6
   plot_opt = 0; % default plot option, false
end

% define x_{1},..., x_{N-1}, total N-1 particles
x = zeros(length(T),N-1); % represent particles x = {x{1},x{2},...,x{N-2},x{N-1}}
x0 = zeros(length(T),1);x_N =zeros(length(T),1); % two boundary particles, are fixed at 0
X = [x0 x x_N]; % position for all N+1 particles
dX = X; ddX = X; % define velocity and acceleration for all particles 

% first N-1 normal modes freq
wk = 2*sin((pi.*(1:5))./(2*N));
Qk = zeros(length(T),length(wk)); % K linear modes
Ek = zeros(length(T),length(wk)); % K modes total energy (kinematic + potential energy)

% initialize
A = zeros(N-1,N-1); B = zeros(N-1,1);
for jj = 1:N-1
A(jj,:) = sqrt(2/N)*sin((pi*jj*(1:N-1))/N);
   if jj == 1
      switch scenario
         case 1 
            B(jj) = 1*sqrt(2)/wk(jj);
         case 2
            B(jj) = sqrt(2*0.75/wk(jj)^2);
      end
   end
    if jj == 3 && scenario == 2
        B(jj) = sqrt(2*0.25/wk(jj)^2);
    end
end

x_0 = A\B; % solve for initial positions
% x_0 = sin(pi*(1:N-1)/(N+1));
X(1,2:N) = x_0';
for kk = 1:length(wk)
   Qk(1,kk) = sqrt(2/N)*sum(X(1,2:end-1).*sin(pi*kk*(1:N-1)/N));
   Ek(1,kk) = 0.5*(((Qk(1,kk) - Qk(1,kk))/dt)^2 + wk(kk)^2*Qk(1,kk)^2);
end

% solving the model
% main time loop
for ii = 2:length(T) 
    
    % main particle loop
    for  jj = 1 : N-1
        if (alpha == 0) && (beta == 0)   % linear model
           ddX(ii, jj+1) = (X(ii-1,jj+2)-2*X(ii-1,jj+1)+X(ii-1,jj));
        elseif alpha~=0     % alpha quadratic nonlinearity model
           ddX(ii, jj+1) = (X(ii-1,jj+2)-2*X(ii-1,jj+1)+X(ii-1,jj)) + alpha *((X(ii-1,jj+2)-X(ii-1,jj+1))^2-(X(ii-1,jj+1)-X(ii-1,jj))^2);       
        elseif beta~=0  % beta cubic nonlinearity model
           ddX(ii, jj+1) = (X(ii-1,jj+2)-2*X(ii-1,jj+1)+X(ii-1,jj)) + beta *((X(ii-1,jj+2)-X(ii-1,jj+1))^3-(X(ii-1,jj+1)-X(ii-1,jj))^3); 
        end
    end
    
    dX(ii,:) = dX(ii-1,:) + ddX(ii,:)*dt;
    X(ii,:) = X(ii-1,:) + dX(ii,:)*dt; 

    % calculate K modes Energy
    for kk = 1:length(wk)
        Qk(ii,kk) = sqrt(2/N)*sum(X(ii,2:end-1).*sin(pi*kk.*(1:N-1)/N));
        Ek(ii,kk) = 0.5*(((Qk(ii,kk) - Qk(ii-1,kk))/dt)^2 + wk(kk)^2*Qk(ii,kk)^2);
    end
    
end


if plot_opt
close;    
figure;
plot(T,Ek(:,1),'r-','LineWidth',1);hold on
plot(T,Ek(:,2),'k-','LineWidth',1);
plot(T,Ek(:,3),'b-','LineWidth',1);
plot(T,Ek(:,4),'m-','LineWidth',1);
plot(T,Ek(:,5),'g-','LineWidth',1);
xlabel('t');ylabel('E_{k}');
legend('mode 1','mode 2','mode 3','mode 4','mode 5')
string = strcat('N=',num2str(N),', alpha=',num2str(alpha),', beta=',num2str(beta),', scenario=',num2str(scenario));
title(string)

end

delta = a*ones(length(T),1);
figure(2)
for i = 1:N+1
    X(:,i) = X(:,i) + i*delta;
end
for j=1:length(T)
    plot(X(j,:),zeros(N+1,1),'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    drawnow;
end

end
