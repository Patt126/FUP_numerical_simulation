function solutionErrors(solutionCord, exactCord,solutionEnergy,nameMethod, dt,t_max,T,MS_cost)

    % Compute error for coordinates
    error_cord = zeros(length(T),6);
    rel_err_cord_inf = zeros(1,6);
    rel_err_cord_ft = zeros(1,6);
    rel_err_cord_l2 = zeros(1,6);

    for i=1:6
        error_cord(:,i) = solutionCord(:,i) - exactCord(:,i)
    
    % Compute relative errors
        rel_err_cord_inf(i) = norm(error_cord(:,i), inf) / norm(exact(:,i), inf);
        rel_err_cord_ft(i) = abs(error_cord(end,i)) / abs(exact(end,i));
        rel_err_cord_l2(i) = norm(error_cord(:,i)) / norm(exact(:,i));
    end
    error_energy = solutionEnergy - 1;
    rel_err_energy_l2 = norm((error_energy))/sqrt(t_max);
    rel_err_energy_ft = abs(error_energy(end));
    rel_err_energy_inf = norm(error_energy,inf);

    % Create or open the CSV file
    file_name = strcat(num2str(MS_cost),'_error.csv');
    
    if exist(file_name, 'file') == 0
        % File does not exist, create it and write headers
        headers = {'Name', 'dt', 'Inf x1', 'L2 x1', 'FT x1','inf x5', 'L2 x5', 'FT x5','Inf x16', 'L2 x16', 'FT x16','Inf Q1', 'L2 Q1', 'FT Q1','inf Q3', 'L2 Q3', 'FT Q3','Inf Q5', 'L2 Q5', 'FT Q5','Inf E', 'L2 E', 'FT E'};
        csvwrite(file_name, headers);
    end

    % Append the results to the CSV file
    result_row = {nameMethod, dt, rel_err_cord_inf(1),rel_err_cord_l2(1), rel_err_cord_ft(1),...
       rel_err_cord_inf(2),rel_err_cord_l2(2), rel_err_cord_ft(2),...
       rel_err_cord_inf(3),rel_err_cord_l2(3), rel_err_cord_ft(3),...
       rel_err_cord_inf(4),rel_err_cord_l2(4), rel_err_cord_ft(4),...
       rel_err_cord_inf(5),rel_err_cord_l2(5), rel_err_cord_ft(5),...
       rel_err_cord_inf(6),rel_err_cord_l2(6), rel_err_cord_ft(6),...
       rel_err_energy_inf,rel_err_energy_l2,rel_err_energy_ft};
    dlmwrite(file_name, result_row, '-append', 'delimiter', ',');

end
