function [rel_err_x_inf, rel_err_energy_inf] = Error_method(x, total_energy, N, alpha, initial_condition , MS_cost, T)

[T_sol, x_exact] = odeSolver(N, alpha, initial_condition, MS_cost, T);

error_x = x - x_exact(:,14);
error_E = total_energy -1;

rel_err_x_inf = norm(error_x, inf);
rel_err_energy_inf = norm(error_E,inf);

end



