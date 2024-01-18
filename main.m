clear all
close all
format long


N = 32;
alpha = 0.1;
n_mode = 5;
t_max = 5*1e3;
number_step = 1e6;
MS_cost = 1; %ratio K/M mass spring constant

[energy_k, x, T_sol] = odeSolver(N, alpha, n_mode, MS_cost, t_max, number_step);

figure(1)
plot(T_sol,energy_k(:,1),'r-','LineWidth',1);hold on
plot(T_sol,energy_k(:,2),'k-','LineWidth',1);
plot(T_sol,energy_k(:,3),'b-','LineWidth',1);
plot(T_sol,energy_k(:,4),'m-','LineWidth',1);
plot(T_sol,energy_k(:,5),'g-','LineWidth',1);
xlabel('t');ylabel('E_{k}');
legend('mode 1','mode 2','mode 3','mode 4','mode 5')

figure(2) 
plot(T_sol,total_energy);
