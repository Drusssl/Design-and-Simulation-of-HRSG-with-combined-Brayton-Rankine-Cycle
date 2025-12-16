% Problem 6: Dynamic Energy-Entropy Cycle Simulation

% Time grid definition
t_end = 3200;
n_points = 1000;
time_vec = linspace(0, t_end, n_points);

%% System Definitions
% Vectorized anonymous functions for temperature and heat
Func_Th = @(t) 900 - 300.*exp(-0.0008.*t);
Func_Tc = @(t) 300 + 40.*sin(0.002.*t);
Func_Qin = @(t) 20000 .* (1 + 0.3.*sin(0.003.*t));

%% Calculation of Profiles
% Evaluate parameters across the entire time vector at once
Th_profile = Func_Th(time_vec);
Tc_profile = Func_Tc(time_vec);
Qin_profile = Func_Qin(time_vec);

% Calculate instantaneous efficiency and power
Efficiency_inst = 1 - (Tc_profile ./ Th_profile);
Power_inst = Efficiency_inst .* Qin_profile;

%% Work Output Integration
% Numerical integration of Power(t) dt
Total_Work = trapz(time_vec, Power_inst);

fprintf('Total Work Output: %.4f kJ\n', Total_Work);

%% Entropy Generation Analysis
% Evaluation using the given relation: S_dot = Q/Th - Q/Tc
S_gen_rate = (Qin_profile ./ Th_profile) - (Qin_profile ./ Tc_profile);

% Visualization
figure(1)
plot(time_vec, S_gen_rate, 'r-', 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Entropy Generation Rate \dot{S}_{gen} (kW/K)')
title('Transient Entropy Generation Profile')
grid on