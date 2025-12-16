% Problem 2: Energy Balance for Closed System

Temp_init = 300;
% Calculate initial internal energy using given u(T) formula
U_init = 450 + 1.1*Temp_init + 0.0012*Temp_init^2;

%% Part 1: Integrate du/dt
time_domain = [0 4000];

% Defined u(T) inversion inline for the ODE function
% Roots of quadratic 0.0012*T^2 + 1.1*T + (450 - u) = 0
% T = (-1.1 + sqrt(1.1^2 - 4*0.0012*(450-u))) / (2*0.0012)
dudt_func = @(time, U_val) 5000*exp(-0.002*time) + ...
    1500*(1 - exp(-0.01 * ((-1.1 + sqrt(1.1^2 - 4*0.0012*(450 - U_val))) / 0.0024)));

[time_vec, U_profile] = ode45(dudt_func, time_domain, U_init);

figure(1)
plot(time_vec, U_profile, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Internal Energy (kJ/kg)')
title('Internal Energy Evolution')
grid on

%% Part 2: Recover Temperature History
% Vectorized calculation (no loop needed)
Temp_history = (-1.1 + sqrt(1.1^2 - 4*0.0012.*(450 - U_profile))) ./ 0.0024;

figure(2)
plot(time_vec, Temp_history, 'r', 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Temperature (K)')
title('Temperature Profile')
grid on

%% Part 3: Crossover Point Analysis
% Calculate heat terms separately to find intersection
Q_ext_vals = 5000 .* exp(-0.002 .* time_vec);
Reaction_vals = 1500 .* (1 - exp(-0.01 .* Temp_history));

% Find index where reaction heat exceeds external heating
idx_cross = find(Reaction_vals >= Q_ext_vals, 1);

if ~isempty(idx_cross)
    crossover_time = time_vec(idx_cross);
    fprintf('Reaction heat surpasses external heating at t approx: %.2f s\n', crossover_time);
else
    fprintf('Crossover point not reached within simulation time.\n');
end