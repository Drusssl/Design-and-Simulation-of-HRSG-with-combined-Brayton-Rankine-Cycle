% Problem 8: Optimal Heating Profile Minimizing Entropy Generation

function Problem8_Optimization()
    % Process Parameters
    Heat_Capacity = 2.5;          
    Temp_Boundary = 300;          
    Temp_Start = 300;          
    Time_End = 2000;         
    Total_Energy_Input = 5e5;        
    
    % Discretization
    Num_Segments = 100;
    Time_Step = Time_End / Num_Segments;
    
    % Initial Guess (Uniform Heating)
    Heat_Rate_Guess = (Total_Energy_Input / Time_End) * ones(Num_Segments, 1);
    
    % Optimization Setup
    % Constraints: No inequality (A,b), Eq constraint for total energy
    LB = zeros(Num_Segments, 1); % Heat rate cannot be negative
    UB = [];
    
    Options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    
    % Solve
    [Heat_Profile_Opt, ~] = fmincon(@(q) Obj_Entropy_Gen(q, Heat_Capacity, Temp_Start, Temp_Boundary, Time_Step), ...
                                    Heat_Rate_Guess, [], [], [], [], LB, UB, ...
                                    @(q) Constr_Total_Heat(q, Total_Energy_Input, Time_Step), Options);
    
    %% Post-Processing
    Heat_Profile_Uniform = Heat_Rate_Guess;
    
    % Reconstruct Temperature Profiles
    % T(t) = T0 + Integral(q/c dt)
    Temp_Profile_Opt = Temp_Start + cumsum(Heat_Profile_Opt) .* (Time_Step / Heat_Capacity);
    Temp_Profile_Uniform = Temp_Start + cumsum(Heat_Profile_Uniform) .* (Time_Step / Heat_Capacity);
    
    % Prepend initial condition for plotting
    Time_Vector = linspace(0, Time_End, Num_Segments);
    
    % Visualization
    figure(1)
    subplot(2,1,1)
    plot(Time_Vector, Heat_Profile_Uniform, 'b--', 'LineWidth', 1.5)
    hold on
    plot(Time_Vector, Heat_Profile_Opt, 'r-', 'LineWidth', 1.5)
    ylabel('Heat Rate q(t) (kW)')
    legend('Uniform Heating', 'Optimal Heating')
    title('Control Profile Comparison')
    grid on
    
    subplot(2,1,2)
    plot(Time_Vector, Temp_Profile_Uniform, 'b--', 'LineWidth', 1.5)
    hold on
    plot(Time_Vector, Temp_Profile_Opt, 'r-', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    title('Temperature Evolution')
    grid on
    
    % Compare Objectives
    S_Gen_Uniform = Obj_Entropy_Gen(Heat_Profile_Uniform, Heat_Capacity, Temp_Start, Temp_Boundary, Time_Step);
    S_Gen_Optimal = Obj_Entropy_Gen(Heat_Profile_Opt, Heat_Capacity, Temp_Start, Temp_Boundary, Time_Step);
    
    fprintf('Entropy Generation (Uniform): %.4f kJ/K\n', S_Gen_Uniform);
    fprintf('Entropy Generation (Optimal): %.4f kJ/K\n', S_Gen_Optimal);
    fprintf('Improvement: %.2f%%\n', (S_Gen_Uniform - S_Gen_Optimal)/S_Gen_Uniform * 100);

end

%% Local Optimization Functions

function S_Gen = Obj_Entropy_Gen(q_vec, C, T_init, T_b, dt)
    % Vectorized temperature calculation
    % current T at step i is approx T_prev + dq
    % Using simple forward accumulation
    Temp_Current = T_init + [0; cumsum(q_vec(1:end-1))] .* (dt / C);
    
    % Entropy Gen Rate = q/T - q/Tb
    % Total S_gen = Sum(Rate * dt)
    Params_Integration = (q_vec ./ Temp_Current) - (q_vec ./ T_b);
    S_Gen = sum(Params_Integration) * dt;
end

function [c_ineq, c_eq] = Constr_Total_Heat(q_vec, Q_target, dt)
    c_ineq = [];
    % Equality: Integral(q dt) = Q_target
    Current_Heat = sum(q_vec) * dt;
    c_eq = Current_Heat - Q_target;
end