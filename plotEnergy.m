function plotEnergy(solution,T,N,n_mode) %aggiungere n_mode 

mode_k = zeros(length(T),N); % modes
energy_k = zeros(length(T),N); %modes' energy
total_energy = zeros(length(T),1);
total_energy(1) = 1;
omega_k = 2*sin(pi.*(1:N)./(2*(N+1))); %frequencies

for i = 1:N
    mode_k(1,i) = sqrt(2/(N+1))*sum(solution(1,N+1:2*N).*sin(pi*i*(1:N)/(N+1)));
    energy_k(1,i) = 0.5*(omega_k(i)*mode_k(1,i))^2;
end

for t = 2:length(T)-1
    for k = 1:N
        mode_k(t,k) = sqrt(2/(N+1))*sum(solution(t,N+1:2*N).*sin(pi*k.*(1:N)/(N+1)));
        energy_k(t,k) = 0.5*((mode_k(t,k)-mode_k(t-1,k))/(T(t)-T(t-1)))^2 + 0.5*(omega_k(k)*mode_k(t,k))^2; 
    end
    total_energy(t) = sum(energy_k(t,:));
end
% Total energy plot
figure();
plot(T,total_energy,'b-','LineWidth',1);
xlabel('time');
ylabel('Total Energy');
% plot of n_modes
figure();
for i = 1:n_mode
    coloreCurva = randi([0, 255], 1, 3);
    coloreCurva = coloreCurva / 255;
    plot(T,energy_k(:,i),'Color',coloreCurva,'LineWidth',1);
    hold on;
end
hold off;
xlabel('time');
ylabel('Energy of modes');
end

