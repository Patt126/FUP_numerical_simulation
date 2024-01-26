function [energy_k,total_energy] = plotEnergy(x, v, T, N, A, MS_cost, n_mode, method_name) 
dt = T(2)-T(1);
mode_k = zeros(length(T),N); 
speed_k = zeros(length(T),N); 
energy_k = zeros(length(T),N); 
total_energy = zeros(length(T),1);
time_average = zeros(length(T),n_mode);
omega_k = 2 * sqrt(MS_cost) * sin(pi .* (1:N)./(2*(N+1))); %frequencies

for t = 1:length(T)
    mode_k(t,:) = (A*x(t,:)')';
    speed_k(t,:) = (A*v(t,:)')';
    energy_k(t,:) = (0.5 * ((speed_k(t,:))).^2 + 0.5 * (omega_k(1,:) .* mode_k(t,:)).^2); 
    total_energy(t) = sum(energy_k(t,:));
      for i = 1:n_mode
          time_average(t,i) = sum(energy_k(1:t,i))/t;
      end
end
%
% plot of n_modes
%
figure(1);
for i = 1:n_mode
    coloreCurva = randi([0, 255], 1, 3);
    coloreCurva = coloreCurva / 255;
    plot(T, energy_k(:,i), 'Color', coloreCurva, 'LineWidth', 1);
    hold on;
end
hold off;
s = strcat('Energy modes  ',method_name);
title(s);
xlabel('time');
ylabel('Energy of modes');
grid on
%
% Show violation of equipartition
%
figure(2);
for i = 1:n_mode
    coloreCurva = randi([0, 255], 1, 3) / 255;
    plot(T, time_average(:,i)/energy_k(1,1), 'Color', coloreCurva, 'LineWidth', 2);
    hold on;
end

% Add horizontal line at -1/32
hline = line([T(1), T(end)], [1/32, 1/32], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2);


hold off;
s = strcat('Temporal average of energy  ',method_name);
title(s);
xlabel('time');
ylabel('Time average of energy');
grid on
legend("mode 1","mode 2","mode 3","mode 4","mode 4")

end


%{
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
s = strcat('Total Energy error  ', method_name);
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
s = strcat('Energy modes  ', method_name);
title(s);
xlabel('time');
ylabel('Energy of modes');
%saveas(2,s)
end
%}




