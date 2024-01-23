function [energy_k,total_energy_error] = plotEnergy(x,v,T,N,A,n_mode, method_name) %aggiungere n_mode 

mode_k = zeros(length(T),N); % modes
speed_k = zeros(length(T),N);
energy_k = zeros(length(T),N); %modes' energy
total_energy_error = zeros(length(T),1);
total_energy_error(1) = 0;
total_total_energy_error = zeros(length(T),1);
total_total_energy_error = 0;
omega_k = 2*sin(pi.*(1:N)./(2*(N+1))); %frequencies

for t = 1:length(T)
    mode_k(t,:) = (A*x(t,:)')';
    speed_k(t,:) = (A*v(t,:)')';
    energy_k(t,:) = 0.5*((speed_k(t,:))).^2 + 0.5*(omega_k(1,:).*mode_k(t,:)).^2; 
    total_energy_error(t) = abs(sum(energy_k(t,:)) -1);
end

for t = 2:length(T)
total_total_energy_error(t) = (total_total_energy_error(t-1) + total_energy_error(t));
end
total_total_energy_error = total_total_energy_error/length(T);
% Choose a window size
window_size = length(T)/1000;

% Calculate the moving average
moving_average = movmean(total_energy_error, window_size);

% Plot the moving average
figure(3)
plot(total_total_energy_error, 'DisplayName', ['Moving Average (Window Size = ', num2str(window_size), ')']);
legend('show');
xlabel('Time');
ylabel('Avarage Value');
s = strcat('Moving Average Analysis',method_name);
title(s);
grid on;
saveas(3,s)
% Total energy plot (CONSIDER PLOT OUT)
figure(1);
plot(T,total_energy_error,'b-','LineWidth',1);
s = strcat('Total Energy error  ',method_name);
title(s);
xlabel('time')
ylabel('Total  error ');
ylim([0,0.2])
%saveas(1,s)
% plot of n_modes
figure(2);
for i = 1:n_mode
    coloreCurva = randi([0, 255], 1, 3);
    coloreCurva = coloreCurva / 255;
    plot(T,energy_k(:,i),'Color',coloreCurva,'LineWidth',1);
    hold on;
end
hold off;
s = strcat('Energy modes  ',method_name);
title(s);
xlabel('time');
ylabel('Energy of modes');
%saveas(2,s)
end

