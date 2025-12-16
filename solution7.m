% Problem 7: Polytropic Piston Compression

% State A Properties
Press_A = 1;   % bar
Temp_A = 300;  % K
Gas_Const = 0.287; % kJ/kg-K

% Process Parameters
Press_B = 10;  % bar
Poly_m = 1.25;

%% 1. Final Temperature Calculation (Analytical)
% Relation: T2/T1 = (P2/P1)^((m-1)/m)
Temp_B = Temp_A * (Press_B / Press_A)^((Poly_m - 1) / Poly_m);

fprintf('Final Temperature (T_B): %.2f K\n', Temp_B);

%% 2. Internal Energy Change
% u(T) = 500 + 0.8T + 1.5e-3*T^2
Func_IntEnergy = @(t) 500 + 0.8.*t + 1.5e-3.*t.^2;

Delta_U = Func_IntEnergy(Temp_B) - Func_IntEnergy(Temp_A);

%% 3. Work and Heat (Analytical)
% Work boundary = R*(T2 - T1) / (1 - m)
Work_Analytical = Gas_Const * (Temp_B - Temp_A) / (1 - Poly_m);

% First Law: Q - W_by = Delta_U  => Q = Delta_U + W_by
% Note: Since compression, Work_Analytical will be negative.
Heat_Transfer = Delta_U + Work_Analytical; 

fprintf('Internal Energy Change: %.4f kJ/kg\n', Delta_U);
fprintf('Heat Transfer: %.4f kJ/kg\n', Heat_Transfer);
fprintf('Analytical Work: %.4f kJ/kg\n', Work_Analytical);

%% 4. Numerical Work Calculation
% Using integral() as per hint
% First, find Specific Volumes at A and B (v = RT/P)
% Note: Pressure in formula must match R units if using standard SI, 
% but here we have ratio P*v^m = C. 
% For Work = Integral(P dv), we need consistent units. 
% 1 bar = 100 kPa. R = 0.287 kJ/kg-K = 287 J/kg-K.
% Let's stick to P in kPa to get Work in kJ.

P_A_kpa = Press_A * 100;
P_B_kpa = Press_B * 100;

Vol_A = (Gas_Const * Temp_A) / P_A_kpa;
Vol_B = (Gas_Const * Temp_B) / P_B_kpa;

% Polytropic Constant K = P * v^m
Poly_K = P_A_kpa * (Vol_A ^ Poly_m);

% Function P(v) = K / v^m
Func_Press_Vol = @(v) Poly_K ./ (v.^Poly_m);

Work_Numerical = integral(Func_Press_Vol, Vol_A, Vol_B);

fprintf('Numerical Work (integral): %.4f kJ/kg\n', Work_Numerical);