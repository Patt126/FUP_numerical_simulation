function [t, u] = mirkf45(N, alpha, initial_condition, MS_cost, T, dt, h0, hmin, TOL)
    facmax = 5;
    fac = 0.9;
    u(:, 1) = initial_condition;
    t(1) = T(1);
    h = h0;
    i = 1;

    m = length(initial_condition);
    orden = 4;  % Order of the method
    c = [0, 1/4, 3/8, 12/13, 1, 1/2];
    A = [0, 0, 0, 0, 0, 0; ...
        1/4, 0, 0, 0, 0, 0; ...
        3/32, 9/32, 0, 0, 0, 0; ...
        1932/2197, -7200/2197, 7296/2197, 0, 0, 0; ...
        439/216, -8, 3680/513, -845/4104, 0, 0; ...
        -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
    b1 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    b2 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];

    K = zeros(m, 6);

    while (t(i) < T(end))
        for j = 1:6
            K(:, j) = feval(@F, t(i) + h * c(j), u(:, i) + h * (K * A(j, :).'), alpha, MS_cost, N);
        end

        % Tentative solution
        phi1 = K * b1.';
        z = u(:, i) + h * phi1;

        % Calculation of the error
        phi2 = K * b2.';
        ERR = norm(phi1 - phi2);

        if (ERR <= TOL)
            t = [t, t(i) + h];  % Use z obtained with the method of lower order.
            u = [u, z];  % Expand the array one by one because we don't know the final size.
            i = i + 1;
        end
           hmax = 0.1;
        h = min([hmax, h * min([facmax, fac * (TOL / ERR) ^ (1 / orden)])]);

        if (h < hmin)
            disp('Error: the step is smaller than hmin.')
            return
        end
    end
end

% Define your F and I functions here
function result = F(t, x, alpha, MS_cost, N)
    e = ones(N, 1);
    linear_part = spdiags([e -2*e e], -1:1, N, N);
    next_quadratic = spdiags([e -e], 0:1, N, N);
    prev_quadratic = spdiags([-e e], -1:0, N, N);

    result = sqrt(MS_cost) * linear_part * x + alpha * ((next_quadratic * x).^2 - (prev_quadratic * x).^2);
end

function result = I(t, state)
    N = length(state);
    result = diag(ones(N, 1)) * state;
end
