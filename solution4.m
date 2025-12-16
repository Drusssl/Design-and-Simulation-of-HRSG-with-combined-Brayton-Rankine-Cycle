% Problem 4: Open System Energy and Entropy Tracking

Temp_in = 310;
Temp_out = 670;
Temp_boundary = 300;

% Thermodynamic Properties
func_Enthalpy = @(t) 300 + 2.5.*t + 0.0007.*t.^2;
func_Entropy = @(t) 2.0.*log(t) + 0.001.*t;

Val_h_in = func_Enthalpy(Temp_in);
Val_h_out = func_Enthalpy(Temp_out);
Val_s_in = func_Entropy(Temp_in);
Val_s_out = func_Entropy(Temp_out);

Delta_h = Val_h_out - Val_h_in;
Delta_s = Val_s_out - Val_s_in;

%% Analysis of Feasible Region
% Range of Heat Transfer Q_dot from 20 to 100 kW
Q_rate_vals = linspace(20, 100, 200);

% 1. Mass flow required by Energy Balance (1st Law): Q = m * dh
M_flow_energy = Q_rate_vals ./ Delta_h;

% 2. Mass flow constraint by Entropy Balance (2nd Law): S_gen >= 0
% m * ds - Q/Tb >= 0  --> m >= Q / (Tb * ds)
M_flow_entropy_limit = Q_rate_vals ./ (Temp_boundary * Delta_s);

% Plotting the comparison
figure(1)
plot(Q_rate_vals, M_flow_energy, 'b-', 'LineWidth', 1.5);
hold on
plot(Q_rate_vals, M_flow_entropy_limit, 'r--', 'LineWidth', 1.5);
xlabel('Heat Rate \dot{Q} (kW)');
ylabel('Mass Flow Rate \dot{m} (kg/s)');
legend('Operating Line (1st Law)', 'Min \dot{m} (2nd Law Limit)', 'Location', 'northwest');
title('Feasible Operating Zone Analysis');
grid on

%% Feasibility Check
% Check if the operating line is ever above the minimum limit
% Since Delta_s/Delta_h is constant, we compare slopes.
Slope_Energy = 1 / Delta_h;
Slope_Entropy = 1 / (Temp_boundary * Delta_s);

if Slope_Energy < Slope_Entropy
    fprintf('Result: No feasible mass flow exists.\n');
    fprintf('Reason: The required mass flow for energy balance is below the minimum required for positive entropy generation.\n');
else
    fprintf('System is feasible.\n');
end

fprintf('\nDiscussion:\n');
fprintf('The enthalpy change fixes the operating slope (First Law).\n');
fprintf('The entropy change sets the lower bound for mass flow (Second Law).\n');
fprintf('Here, the energy line lies in the forbidden region (S_gen < 0).\n');