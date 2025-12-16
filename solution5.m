% Problem 5: Exergy Analysis with Temperature-Dependent Properties

Temp_in = 350;
Temp_out = 900;
Temp_dead = 298;

%% Part 1: Numerical Entropy Change
% Cp(T) = 1200 + 0.4T - 1.2e-4*T^2 (J/kg-K)
% Entropy differential: ds = Cp(T)/T dT
Func_Entropyter = @(t) (1200 + 0.4.*t - 1.2e-4.*t.^2) ./ t;

% Using adaptive quadrature (integral) as per hint
Delta_S_sys = integral(Func_Entropyter, Temp_in, Temp_out);

fprintf('Entropy Change of Stream: %.2f J/kg-K\n', Delta_S_sys);

%% Part 2: Exergy Destruction Scenarios
% The problem asks to estimate exergy destruction for specific irreversibility levels
% Assuming "Irreversibility of X%" implies S_gen = X% of Delta_S
Irr_levels = [0.02, 0.10]; 

fprintf('\nExergy Destruction Analysis:\n');
for k = 1:length(Irr_levels)
    S_gen_k = Irr_levels(k) * Delta_S_sys;
    X_dest_k = Temp_dead * S_gen_k; % Gouy-Stodola Theorem
    fprintf('At %.0f%% Irreversibility: %.4f kJ/kg\n', Irr_levels(k)*100, X_dest_k/1000);
end

%% Part 3: Sensitivity Analysis
Fraction_Irr = linspace(0.02, 0.10, 100);
X_dest_vector = Temp_dead .* (Fraction_Irr .* Delta_S_sys);

figure(1)
plot(Fraction_Irr * 100, X_dest_vector / 1000, 'b-', 'LineWidth', 1.5)
xlabel('Heater Irreversibility (%)')
ylabel('Exergy Destruction (kJ/kg)')
title('Impact of Irreversibility on Exergy Destruction')
grid on