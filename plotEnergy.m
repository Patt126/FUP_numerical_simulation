function [energy_k,total_energy] = plotEnergy(x, v, T, N, A, MS_cost, n_mode, method_name) 
dt = T(2)-T(1);
mode_k = zeros(length(T),N); 
speed_k = zeros(length(T),N); 
energy_k = zeros(length(T),N); 
total_energy = zeros(length(T),1);
time_average = zeros(length(T),n_mode);
omega_k = 2 * sqrt(MS_cost) * sin(pi .* (1:N)./(2*(N+1))); %frequencies

for t = 1:length(T)
    mode_k(t,:) = (A * x(t,:)')';
    speed_k(t,:) = (A * v(t,:)')';
    energy_k(t,:) = 0.5 * ((speed_k(t,:))).^2 + 0.5 * (omega_k(1,:) .* mode_k(t,:)).^2; 
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
s = strcat('Energy modes  ', method_name);
title(s);
xlabel('time');
ylabel('Energy of modes');
%
% Show violation of equipartition
%
figure(2);
for i = 1:n_mode
    coloreCurva = randi([0, 255], 1, 3);
    coloreCurva = coloreCurva / 255;
    plot(T, time_average(:,i), 'Color', coloreCurva, 'LineWidth', 1);
    hold on;
end
hold off;
s = strcat('Temporal average of energy  ', method_name);
title(s);
xlabel('time');
ylabel('Time average of energy');

end



