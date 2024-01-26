clear all
close all
format long
%
file_path = '1_error.csv';

% Use readtable
table_data = readtable(file_path);
% Per accedere a una colonna table_data.nomeColonna
% Ricordare che siamo andati da dt max a dt min quindi bisogna sempre
% invertire gli array con flip()!

number_of_methods = 5;
number_of_dt = 7;

dt = flip(table_data.dt);
dt = dt(1:number_of_dt);

figure(1)
title("Linf X16");
hold on 
grid on
start_index = 1;
end_index = number_of_dt;
data = (table_data.InfX16); %put here the data to plot

for i = 1:number_of_methods
    data_current = flip(data(start_index:end_index)); %reverse data order as dt
    plot(dt,data_current, 'LineWidth',2,'Marker', 'o') %plot different data
    start_index = start_index + number_of_dt ;
    end_index = end_index + number_of_dt ;
end
hold off
legend("Euler","Leapfrog","Nystrom","ruth","rk4");

%evaluate peff for each method on different coordinates 
start_index = 22;
end_index = 28;
energy_inf = table_data.InfE;
energy_inf = flip(energy_inf(start_index:end_index));
energy_l2 = table_data.L2E;
energy_l2 = flip(energy_l2(start_index:end_index));
Q1_inf = table_data.infQ3;
Q1_inf = flip(Q1_inf(start_index:end_index));
Q1_l2 = table_data.L2Q3;
Q1_l2 = flip(Q1_l2(start_index:end_index));

peff_E_inf = zeros(1,number_of_dt);
peff_E_l2 = zeros(1,number_of_dt);
peff_Q_inf = zeros(1,number_of_dt);
peff_Q_l2 = zeros(1,number_of_dt);


for i = 2:number_of_dt
    peff_E_inf(i) = log2(energy_inf(i)/energy_inf(i-1));
    peff_E_l2(i) = log2(energy_l2(i)/(sqrt(2) * energy_l2(i-1)));
    peff_Q_l2(i) = log2(Q1_l2(i)/(sqrt(2) * Q1_l2(i-1)));
    peff_Q_inf(i) = log2(Q1_inf(i)/(sqrt(2) * Q1_inf(i-1)));
    
end

