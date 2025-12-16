% Problem 10: Entropy-Constrained Design of a Heat Pump Stage

function Problem10_HeatPumpDesign()
    % System Parameters
    Temp_Cold = 270; % K
    Temp_Hot  = 320; % K
    Alpha_Loss = 0.02;
    K_Work = 50;

    % Optimization Bounds for Pressure Ratio (rp)
    Rp_Lower = 1.1;
    Rp_Upper = 6.0;
    Rp_Guess = 2.0;

    % Optimization Setup
    % Objective: Maximize COP (equivalent to Minimizing -COP)
    Obj_Func = @(rp) -Func_COP(rp, Temp_Cold, Temp_Hot, Alpha_Loss);
    
    % Constraint: Entropy Generation <= 0.05
    Constraint_Func = @(rp) NonLinear_Constraints(rp, Temp_Cold, Temp_Hot, Alpha_Loss, K_Work);
    
    Options = optimoptions('fmincon', 'Display', 'notify', 'Algorithm', 'sqp');

    % Solve
    [Rp_Optimal, Neg_COP_Opt] = fmincon(Obj_Func, Rp_Guess, [], [], [], [], ...
                                        Rp_Lower, Rp_Upper, Constraint_Func, Options);

    % Extract Results at Optimum
    COP_Max = -Neg_COP_Opt;
    S_Gen_Opt = Func_EntropyGen(Rp_Optimal, Temp_Cold, Temp_Hot, Alpha_Loss, K_Work);

    fprintf('--- Optimization Results ---\n');
    fprintf('Optimal Pressure Ratio (rp): %.4f\n', Rp_Optimal);
    fprintf('Maximum COP: %.4f\n', COP_Max);
    fprintf('Entropy Generation at Opt: %.4f kJ/K\n', S_Gen_Opt);

    %% Visualization
    Rp_Vector = linspace(1.1, 6, 500);
    
    % Evaluate profiles
    COP_Profile = Func_COP(Rp_Vector, Temp_Cold, Temp_Hot, Alpha_Loss);
    S_Gen_Profile = Func_EntropyGen(Rp_Vector, Temp_Cold, Temp_Hot, Alpha_Loss, K_Work);
    
    % Identify Feasible Region (S_gen <= 0.05)
    Feasible_Limit = 0.05;
    Is_Feasible = S_Gen_Profile <= Feasible_Limit;

    figure(1)
    
    % Subplot 1: COP Behavior
    subplot(3,1,1)
    plot(Rp_Vector, COP_Profile, 'b-', 'LineWidth', 1.5)
    ylabel('COP')
    title('Performance Characteristics vs Pressure Ratio')
    grid on
    
    % Subplot 2: Entropy Generation Constraint
    subplot(3,1,2)
    plot(Rp_Vector, S_Gen_Profile, 'r-', 'LineWidth', 1.5)
    yline(Feasible_Limit, 'k--', 'Constraint (0.05)', 'LineWidth', 1.2)
    ylabel('S_{gen} (kJ/K)')
    grid on
    
    % Subplot 3: Feasible Operating Region
    subplot(3,1,3)
    plot(Rp_Vector, COP_Profile, 'k:', 'LineWidth', 1) % Full curve phantom
    hold on
    plot(Rp_Vector(Is_Feasible), COP_Profile(Is_Feasible), 'g-', 'LineWidth', 2.5) % Feasible part
    plot(Rp_Optimal, COP_Max, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r') % Opt point
    xlabel('Pressure Ratio r_p')
    ylabel('Feasible COP')
    legend('Infeasible Region', 'Feasible Region', 'Optimal Point', 'Location', 'best')
    grid on

    %% Discussion
    fprintf('\n--- Design Trade-offs ---\n');
    fprintf('1. Higher pressure ratios generally decrease COP due to the loss term alpha*(rp-1)^2/rp.\n');
    fprintf('2. However, compressor work increases with rp. Entropy generation depends on both Heat and Work.\n');
    fprintf('3. The constraint forces us to operate away from the global unconstrained maximum if S_gen is too high.\n');
    fprintf('   (In this specific case, the algorithm finds the best balance within the limit.)\n');
end

%% Local Helper Functions

function val_COP = Func_COP(rp, Tc, Th, alpha)
    % Model provided in assignment [cite: 121]
    % Note: Formula uses Tc in numerator (Refrigeration-like), but we follow the PDF strictly.
    Carnot_Factor = Tc / (Th - Tc);
    Loss_Factor = 1 - alpha .* ((rp - 1).^2) ./ rp;
    val_COP = Carnot_Factor .* Loss_Factor;
end

function val_Work = Func_CompressorWork(rp, k)
    % Wc = k * (sqrt(rp) - 1) [cite: 124]
    val_Work = k .* (sqrt(rp) - 1);
end

function S_gen = Func_EntropyGen(rp, Tc, Th, alpha, k)
    % Calculate Flows
    Current_COP = Func_COP(rp, Tc, Th, alpha);
    Work_In = Func_CompressorWork(rp, k);
    
    % Heat definitions based on COP = Qh / W (for Heat Pump) or Qc / W?
    % The problem asks for "Heat delivered per unit work" which implies Heating COP.
    % Assuming Standard Thermodynamics: W + Qc = Qh
    % If COP defined as Q_delivered / W:
    Q_Delivered = Current_COP .* Work_In; % Qh
    Q_Source = Q_Delivered - Work_In;     % Qc
    
    % Entropy Balance: S_gen = (Q_out/T_out) - (Q_in/T_in)
    % S_gen = Qh/Th - Qc/Tc
    S_gen = (Q_Delivered ./ Th) - (Q_Source ./ Tc);
end

function [c_ineq, c_eq] = NonLinear_Constraints(rp, Tc, Th, alpha, k)
    % Constraint: S_gen <= 0.05
    Current_S_gen = Func_EntropyGen(rp, Tc, Th, alpha, k);
    
    c_ineq = Current_S_gen - 0.05;
    c_eq = [];
end