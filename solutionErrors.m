function [rel_err_inf, rel_err_l2, rel_err_final] = solutionErrors(solution, exact,nameCoord ,nameMethod, dt)

    % Compute error
    error = solution - exact;

    % Compute relative errors
    rel_err_inf = norm(error, inf) / norm(exact, inf);
    rel_err_final = abs(error(end)) / abs(exact(end));
    rel_err_l2 = norm(error) / norm(exact);

    % Create or open the CSV file
    file_name = strcat(nameCoord,'_results.csv');
    
    if exist(file_name, 'file') == 0
        % File does not exist, create it and write headers
        headers = {'Name', 'dt', 'RelativeErrorInf', 'RelativeErrorL2', 'RelativeErrorFinal'};
        csvwrite(file_name, headers);
    end

    % Append the results to the CSV file
    result_row = {name, dt, rel_err_inf, rel_err_l2, rel_err_final};
    dlmwrite(file_name, result_row, '-append', 'delimiter', ',');

end
