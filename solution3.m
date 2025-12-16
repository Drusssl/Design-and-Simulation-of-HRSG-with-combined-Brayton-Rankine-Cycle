% Problem 3: Entropy in Compression with Real-Gas Effects

% Constants
Coeff_const = 1;
Coeff_P = 0.0008;
Coeff_T_inv = -120;
Poly_n = 1.28;
Cp_mass = 1.05; % kJ/kg.K
R_gas = 0.287;  % kJ/kg.K
Temp_inlet = 300;

Press_start = 1;
Press_end = 20;

%% 1. Solve for Temperature Profile T(P)
% Z factor function
Func_Z = @(P, T) Coeff_const + Coeff_P.*P + Coeff_T_inv./T;

% Differential equation dT/dP derived from Polytropic relation
% Derived from (1-n)lnP + n*lnZ + n*lnT = Const
dTdP_eqn = @(P, T) (Func_Z(P,T) .* T.^2 .* ( (Poly_n - 1)./P - Poly_n*Coeff_P./Func_Z(P,T) )) ./ ...
                   (Poly_n .* (Func_Z(P,T).*T - Coeff_T_inv));

% Define a high-resolution pressure grid for integration
Res_points = 200;
Press_grid = linspace(Press_start, Press_end, Res_points);

[Press_sol, Temp_sol] = ode45(dTdP_eqn, Press_grid, Temp_inlet);

figure(1)
plot(Press_sol, Temp_sol, 'b-', 'LineWidth', 1.5)
xlabel('Pressure (bar)')
ylabel('Temperature (K)')
title('Temperature vs Pressure (Polytropic Real Gas)')
grid on

%% 2. Compute Real Gas Entropy Change
% Calculate Z values along the path
Z_values = Func_Z(Press_sol, Temp_sol);

% Entropy change: ds = Cp*dT/T - R*Z*dP/P
% Integrating: Integral(Cp/T dT) - Integral(R*Z/P dP)
Term1_Real = trapz(Temp_sol, Cp_mass ./ Temp_sol);
Term2_Real = trapz(Press_sol, R_gas .* Z_values ./ Press_sol);

Delta_S_Real = Term1_Real - Term2_Real;

fprintf('Real Gas Entropy Change: %.6f kJ/kg-K\n', Delta_S_Real);

%% 3. Compute Ideal Gas Entropy Change
% Ideal case: ds = Cp*dT/T - R*dP/P (Z=1)
Term1_Ideal = trapz(Temp_sol, Cp_mass ./ Temp_sol);
Term2_Ideal = trapz(Press_sol, R_gas ./ Press_sol);

Delta_S_Ideal = Term1_Ideal - Term2_Ideal;

fprintf('Ideal Gas Entropy Change: %.6f kJ/kg-K\n', Delta_S_Ideal);

%% 4. Deviation Analysis
Deviation_pct = abs((Delta_S_Ideal - Delta_S_Real) / Delta_S_Ideal) * 100;
fprintf('Percent Deviation: %.4f%%\n', Deviation_pct);