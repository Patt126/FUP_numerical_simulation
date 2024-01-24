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
data = (table_data.InfX16);

for i = 1:number_of_methods
    data_current = flip(data(start_index:end_index));
    plot(dt,data_current, 'LineWidth',2,'Marker', 'o')
    start_index = start_index + number_of_dt ;
    end_index = end_index + number_of_dt ;
end
hold off
legend("Euler","Leapfrog","Nystrom","ruth","rk4");
