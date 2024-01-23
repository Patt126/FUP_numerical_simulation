function [ x_values, v_values, t_values] = adaptive(N, alpha, initial_condition, MS_cost, T, tol)
    x_values = zeros(length(T), N);
    x_values(1, :) = initial_condition(N+1:2*N);
    
    v_values = zeros(length(T), N);
    v_values(1, :) = initial_condition(1:N);
    
    t_values = zeros(length(T), 1);
    t_values(1) = T(1);
    
    dt_values = zeros(length(T), 1);

    e = ones(N,1);
    linear_part = spdiags([e -2*e e],-1:1,N,N);
    next_quadratic = spdiags([e -e], 0:1 ,N,N);
    prev_quadratic = spdiags([-e e], -1:0 ,N,N);
    
    F = @(t, x) sqrt(MS_cost) * linear_part * x + alpha * ((next_quadratic * x).^2 - (prev_quadratic * x).^2);
    I = @(t, state) diag(ones(N, 1)) * state;

    idx = 1;
    while t_values(idx) < T(end)
        % Perform a single step using DP method
        [t_next, x_next, v_next, dt_next] = dp_step(t_values(idx), x_values(idx, :)', v_values(idx, :)', F, I, MS_cost, alpha, tol, N);

        % Update the arrays
        idx = idx + 1;
        t_values(idx) = t_next;
        x_values(idx, :) = x_next';
        v_values(idx, :) = v_next';
        dt_values(idx) = dt_next;
    end
end

function [t_next, x_next, v_next, dt_next] = dp_step(t, x, v, F, I, MS_cost, alpha, tol, N)
    % Dormand-Prince parameters
    c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
    a = [
        0, 0, 0, 0, 0, 0, 0;
        1/5, 0, 0, 0, 0, 0, 0;
        3/40, 9/40, 0, 0, 0, 0, 0;
        44/45, -56/15, 32/9, 0, 0, 0, 0;
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0;
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0
    ];
    b = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
    b_hat = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];

    dt = 0.01; % Initial step size
    max_dt = 0.1; % Maximum step size

    k = zeros(N, length(c));

    while true
        % Compute the slopes using DP coefficients
        for i = 1:length(c)
           k(:, i) = dt * F(t + c(i)*dt, x + k * a(i, :)')';

        end

        % Compute the next state using weighted sums
        x_next = x + k * b';
        v_next = v + I(t, k(:, 1));

        % Compute the estimated error
        error = norm(b_hat' * k', 2);


        % Update the time and check if the error is within tolerance
        t_next = t + dt;
         
        if error < tol || dt <= tol
            break; % Accept the step
        else
            % Adjust the step size
            dt = min(0.9 * dt * (tol / error)^(1/5), max_dt);
            dt_next = dt;
        end
    end
end
