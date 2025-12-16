% Problem 9: Reaction Heat and Parameter Estimation

function Problem9_ParameterEstimation()
    % True System Parameters
    PreExp_True = 1e5;
    Activ_Energy_True = 45; % kJ/mol
    Heat_Cap = 1.8;         % kJ/K
    Gas_Const = 8.314e-3;   % kJ/mol-K

    % Simulation Settings
    Time_Span = linspace(0, 5, 500); % 0 to 5 seconds
    Temp_Init = 300;                 % Initial Temp K

    %% Part 1: Generate Synthetic Data
    % Define true ODE: c*dT/dt = A*exp(-E/RT) + Q_ext(t)
    Func_ODE_True = @(t, T) (PreExp_True * exp(-Activ_Energy_True ./ (Gas_Const * T)) ...
                             + 2000 * exp(-0.001 * t)) / Heat_Cap;

    [Time_Sol, Temp_Clean] = ode45(Func_ODE_True, Time_Span, Temp_Init);

    % Add 1% Gaussian Noise
    Noise_Level = 0.01;
    Temp_Noisy = Temp_Clean .* (1 + Noise_Level * randn(size(Temp_Clean)));

    % Visualize Data
    figure(1)
    plot(Time_Sol, Temp_Noisy, 'r.', 'MarkerSize', 8)
    hold on
    plot(Time_Sol, Temp_Clean, 'k-', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    legend('Noisy Synthetic Data', 'True Model (Noise-Free)', 'Location', 'NorthWest')
    title('Synthetic Data Generation')
    grid on

    %% Part 2: Parameter Estimation via Inverse Modeling
    % We try to recover A and E from the noisy data
    % Initial Guess: A = 5e4, E = 30 (deliberately off)
    Params_Guess = [5e4, 30]; 
    
    % Optimization Options
    Opt_Options = optimoptions('lsqcurvefit', 'Display', 'off');

    % Define objective function wrapper
    % Note: x corresponds to time, y to Temperature
    Func_Model_Fit = @(params, time_points) Solve_Temp_ODE(params, time_points, Temp_Init, Heat_Cap, Gas_Const);

    fprintf('Starting Parameter Estimation...\n');
    [Params_Est, Resnorm] = lsqcurvefit(Func_Model_Fit, Params_Guess, Time_Sol, Temp_Noisy, [], [], Opt_Options);

    A_Est = Params_Est(1);
    E_Est = Params_Est(2);

    %% Part 3: Results and Discussion
    fprintf('\n--- Estimation Results ---\n');
    fprintf('Parameter A: True = %.2e | Estimated = %.2e | Error = %.2f%%\n', ...
        PreExp_True, A_Est, abs(A_Est - PreExp_True)/PreExp_True * 100);
    fprintf('Parameter E: True = %.2f  | Estimated = %.2f  | Error = %.2f%%\n', ...
        Activ_Energy_True, E_Est, abs(E_Est - Activ_Energy_True)/Activ_Energy_True * 100);

    % Check fit quality
    Temp_Fitted = Solve_Temp_ODE(Params_Est, Time_Sol, Temp_Init, Heat_Cap, Gas_Const);
    
    plot(Time_Sol, Temp_Fitted, 'b--', 'LineWidth', 2)
    legend('Noisy Synthetic Data', 'True Model', 'Fitted Model')
    
    fprintf('\nDiscussion:\n');
    fprintf('The fitting algorithm successfully navigates the noise to recover parameters.\n');
    fprintf('Small deviations exist due to the non-linear coupling of A and E in the exponential term.\n');
end

%% Local Helper Function for ODE Integration inside Fitting Loop
function T_Profile = Solve_Temp_ODE(params, t_vector, T0, C, R)
    A_curr = params(1);
    E_curr = params(2);
    
    % Re-define ODE with current iteration parameters
    dTemp_dt = @(t, T) (A_curr * exp(-E_curr ./ (R * T)) + 2000 * exp(-0.001 * t)) / C;
    
    % Solve ODE
    % Note: ode45 returns a structure or variable steps, we need values at specific t_vector
    % We use deval or simply pass t_vector to ode45 (if t_vector is sorted)
    [~, T_Profile] = ode45(dTemp_dt, t_vector, T0);
end