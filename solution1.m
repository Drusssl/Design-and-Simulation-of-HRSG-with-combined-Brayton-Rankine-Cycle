% Problem 1: Internal Energy Estimation

Temp_start = 320;
Temp_end = 820;

%% Symbolic Solution
syms Temp_sym
Cv_eqn = 700 + 0.35*Temp_sym - 2e-4*Temp_sym^2;
Delta_U_exact = int(Cv_eqn, Temp_sym, Temp_start, Temp_end);
Val_Exact = double(Delta_U_exact);

fprintf('Analytical Result (Symbolic): %.6f kJ/kg\n', Val_Exact/1000);

%% Numerical Method 1: integral()
func_Cv = @(t) 700 + 0.35.*t - 2e-4.*t.^2;
Delta_U_integral = integral(func_Cv, Temp_start, Temp_end);

fprintf('Numerical Result (integral): %.6f kJ/kg\n', Delta_U_integral/1000);

%% Numerical Method 2: Grid Convergence with trapz()
n_grid = 2;
tolerance = 0.001; 

while true
    T_range = linspace(Temp_start, Temp_end, n_grid);
    Cv_vals = func_Cv(T_range);
    
    Delta_U_approx = trapz(T_range, Cv_vals);
    
    curr_error = abs((Delta_U_approx - Val_Exact) / Val_Exact);
    
    if curr_error <= tolerance
        break;
    end
    n_grid = n_grid + 1;
end

fprintf('Numerical Result (trapz): %.6f kJ/kg\n', Delta_U_approx/1000);
fprintf('Minimum grid points required: %d\n', n_grid);