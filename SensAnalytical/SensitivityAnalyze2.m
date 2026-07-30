function SensitivityAnalyze2(PID)
% Complete cardiovascular sensitivity analysis workflow adapted from software reference framework
   %% 1. Run the core simulation dynamically for the requested patient
   tic
    fprintf('Running analytical LSA for Patient %d...\n', PID);
    cd ..
    
    [~, S_history] = DriverBasic(PID, 0, 0, 1, [0 0 1], 12);
    cd SensAnalytical
    %SolveModel(pars,data);
    %% 2. Retrieve variables from the base workspace
    if ~evalin('base', 'exist(''S_history'', ''var'')')
        error('Could not find S_history in workspace. Ensure additions are saved in DriverBasic.m')
    end

    S_history = evalin('base', 'S_history'); % Dimension: [6 states x 16 parameters x N time steps]
    
    S_time_raw = evalin('base', 'S_time'); % Vector of length N
    
    [num_states, num_params, num_steps] = size(S_history); 
        
    % Directly grab the unlogged raw parameter vector exported from the driver
    cd ../Core
    pars_vector = evalin('base', 'pars_nominal');
    
    % Reconstruct the data structure and logged parameter that senseq
    % expects
    % This matches DriverBasic_sens structure loading workflow
    data_struct = Patient(PID, 0);
    [pars_from_load, parsG, times, Init] = load_global(PID, data_struct);
    data_struct.parsG = parsG;
    data_struct.Init  = Init;
    data_struct.dt    = times.dt;
    data_struct.T     = times.T;
    data_struct.CC    = 5;
    data_struct.computeAnalyticalSens = true; %%%%% HERE BRO %%%%%

    % Reconstruct the logged parameters that senseq expects
    pars_logged = log(pars_vector);

    % % Fetch the patient table from the base workspace memory to get labels
    % patientTable = evalin('base', 'patientTable');

    %% 3. Slice Time to Match S_history Steps Exactly
    len_time = length(S_time_raw);
    
    if len_time == num_steps
        % Perfect match
        S_time = S_time_raw;
    elseif num_steps == len_time + 1
        % S_history has exactly one extra frame at the end (off-by-one loop artifact)
        % Safely truncate the extra history frame to perfectly match the time vector
        S_history = S_history(:, :, 1:len_time);
        %S_scaled = S_scaled(:, :, 1:len_time); % Re-align dimensions for later loops
        num_steps = len_time; % Update tracking variable
        S_time = S_time_raw;
    elseif len_time > num_steps
        % Time vector is longer than history frames; slice from the end safely
        S_time = S_time_raw(end-num_steps+1:end);
    else
        % Fallback: If time vector is significantly shorter than history frames
        if len_time > 1
            dt = mean(diff(S_time_raw));
        else
            dt = 0.01; % Safe default step fallback size (10 ms)
        end
        % Reconstruct timeline matching the step frequency
        S_time = 0:dt:(num_steps-1)*dt;
    end
    
    %% 4. Perform Parameter Relative Scaling (Individually Standardized per State Compartment)
    % Allocate empty tracking array
    S_scaled = zeros(num_states, num_params, num_steps);
    
    for s_idx = 1:num_states
        % Find the absolute peak value achieved by ANY parameter in this state across all timesteps
        state_slice = S_history(s_idx, :, :);
        max_state_val = max(abs(state_slice(:)));
        
        if max_state_val == 0
            max_state_val = 1; % Prevent divide-by-zero errors
        end
        
        for p_idx = 1:num_params
            raw_curve = S_history(s_idx, p_idx, :);
            % Normalize every curve by its specific compartment maximum scalar
            S_scaled(s_idx, p_idx, :) = raw_curve / max_state_val;
        end
    end
    
    % Allocate a 2D matrix to store a score for every state-parameter pair
    Rsens_matrix = zeros(num_states, num_params);
    
    for s_idx = 1:num_states
        for p_idx = 1:num_params
            % Pull out 1D vector of timesteps for this state and this parameter
            time_vector_slice = squeeze(S_scaled(s_idx, p_idx, :));
            % Calculate the regular Euclidean norm over time
            Rsens_matrix(s_idx, p_idx) = norm(time_vector_slice, 2);
        end 
    end     
    
    % Define structural string labels for indexing mappings
    param_labels = {'Rpul', 'Rmval', 'Raval', 'Rsys', 'Rtval', 'Rpval', ...
                    'Cpv', 'Csa', 'Csv', 'Cpa', 'ElvM', 'Elvm', 'ErvM', 'Ervm', 'TC', 'TR'};
    state_names  = {'Pulmonary Veins (Vpv)', 'Left Ventricle (Vlv)', 'Systemic Arteries (Vsa)', ...
                    'Systemic Veins (Vsv)', 'Right Ventricle (Vrv)', 'Pulmonary Arteries (Vpa)'};
    
    %% 5. 3x2 Tiled Ranked Sensitivity Plots
    fig1 = figure('Name', sprintf('Patient %d: State-Specific Ranked Parameter Importance', PID));
    
    % Initialize tiled layout
    t_rank = tiledlayout(3, 2, 'TileSpacing', 'Loose', 'Padding', 'Compact');

    for s_idx = 1:num_states
        nexttile;
        state_vals = Rsens_matrix(s_idx, :);

        % Normalize the importance scores for this specific state to fall
        % between 0 and 1
        max_local_score = max(state_vals);
        if max_local_score == 0, max_local_score = 1; end
        state_vals_normalized = state_vals / max_local_score;

        % Sort just this state's values in descending order
        [sorted_vals, Isens_local] = sort(state_vals_normalized, 'descend');
    
        % Plot points as explicit blue X markers 
        plot(1:16, sorted_vals, 'x', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', '#0072BD');
        
        % Explicitly clip and lock the horizontal X-axis limits to 1-16 to avoid data stretching
        xlim([0.5, 16.5]);
        xticks(1:16);
        ylim([0, 1.05])
        
        % Maps your sorted human-readable labels to the x-axis ticks
        xticklabels(param_labels(Isens_local)); 
        xtickangle(45);
        
        % Plot labels and aesthetic configurations
        ylabel('Rel. Sens. Norm [0-1]', 'FontSize', 8, 'FontWeight', 'bold');
        title(sprintf('%s', state_names{s_idx}), 'FontSize', 10);
        grid on;

        % Only append horizontal X-axis labels on the bottom two cells to maximize layout space
        if s_idx >= 5
            xlabel('Parameters (Ranked)', 'FontSize', 10, 'FontWeight', 'bold');
        end

        
    end 

    % Add a master title spanning the top of tiled layout
    title(t_rank, sprintf('Patient %d: Parameter Importance for Volumetric States', PID), ...
          'FontSize', 12, 'FontWeight', 'Bold')

    %% 6. 3x2 Tiled Time-Series Sensitivity Plots
    % Zero out the time vector so visualization plots explicitly start at 0.0 seconds
    S_time_zeroed = S_time - S_time(1);

    % Open one singular background figure window to contain all 6 plots
    fig_time = figure('Name', sprintf('Patient %d: Time-Varying Relative Sensitivities', PID));
    
    % Initialize a 3-row, 2-column tiled canvas layout
    t = tiledlayout(3, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    
    % Universal line configurations matching your style framework
    line_styles = {'-', '--', ':', '-.'};
    
    for s_idx = 1:num_states
        % Move to the next sub-plot tile coordinate in the 3x2 grid matrix
        nexttile;
        hold on;
        
        % Loop through and draw the curves over time for all 16 parameters
        for p_idx = 1:num_params
            time_varying_curve = S_scaled(s_idx, p_idx, :);
            time_varying_curve = time_varying_curve(:);
            
            % Cycle lines dynamically through your 4 native style types
            style_idx = mod(p_idx - 1, 4) + 1;
            
            % Plot standard linear scale values
            plot(S_time_zeroed, abs(time_varying_curve), ...
                 'LineStyle', line_styles{style_idx}, 'LineWidth', 1);
        end
        
        hold off;
        grid on;
        set(gca, 'FontSize', 12);

        % Force the y-axis scale boundary limits
        ylim([0, 1.05])
        
        % Individual tile descriptions
        ylabel('Rel. Sensitivity', 'FontSize', 10, 'FontWeight', 'bold');
        title(sprintf('%s', state_names{s_idx}), 'FontSize', 10);
        
        % Only append horizontal X-axis labels on the bottom two cells to maximize layout space
        if s_idx >= 5
            xlabel('Time (s)', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end 
    
    % Create one clean master legend column on the far right side of the layout canvas
    lgd = legend(param_labels, 'Orientation', 'vertical', 'FontSize', 12);
    lgd.Layout.Tile = 'east'; 
    
    % Add an overarching header title bar spanning across the top of the dashboard image
    title(t, sprintf('Patient %d: Continuous Relative Parameter Sensitivities Across the Cardiac Cycle', PID), ...
          'FontSize', 12, 'FontWeight', 'Bold');

    %% 7. Archive Data Out to File
    cd ../SensAnalytical/
    
    save('CardioSensResults.mat', 'S_scaled', 'Rsens_matrix');
    fprintf('Data profiles successfully archived to file "CardioSensResults.mat"\n');

    %% 8. Create Overarching Single Parameter Ranking Plot (Using Mette's Code Framework)
    % Execute senseq using finite difference logic
    cd ../Core/
    sens_other = senseq(pars_logged, data_struct);

    [~, N_m] = size(sens_other);
    sens_norm_other = zeros(1, N_m);
    for i = 1:N_m
        sens_norm_other(i) = norm(sens_other(:, i), 2);
    end    
    
    % Normalize relative to the single largest overall parameter rank
    sens_norm_scaled = sens_norm_other / max(sens_norm_other);

    % Sort from highest sensitivity to lowest
    [Rsens, Isens] = sort(sens_norm_scaled, 'descend');

    % Create the single comparison plot figure
    fig2 = figure('Name', sprintf('Patient %d: Analytical Parameter Sensitivity Ranking', PID));
    
    %% Compute overall analytical sensitivities (1x16)

    analytical_overall = zeros(1, num_params);
    
    for p = 1:num_params
        analytical_overall(p) = norm(Rsens_matrix(:,p),2);
    end
    
    % Normalize to match finite-difference scaling
    analytical_overall = analytical_overall ./ max(analytical_overall);
    
    %% Compare the two methods
    
    abs_error = abs(analytical_overall - sens_norm_scaled);
    
    rel_error = 100*abs_error ./ max(abs(sens_norm_scaled),eps);
        
    % Plot points as explicit blue X markers matching the style of the tiled layout
    plot(1:16, Rsens, 'x', 'LineWidth', 4, 'MarkerSize', 20, 'Color', '#0072BD');
    grid on;
    
    % Explicitly clip and lock the horizontal X-axis limits to 1-16 to avoid data stretching
    xlim([0.5, 16.5]);
    xticks(1:16);
    ylim([0, 1.05])
    
    set(gca, 'XTick', 1:16, 'XTickLabel', param_labels(Isens), 'FontSize', 20);
    xtickangle(45);
    
    ylabel('Rel. Volumetric Sensitivity Norm [0-1]', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Parameters (Ranked)', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Patient %d: Trajectory Sensitivity Rankings (Analytical Chain-Rule Method)', PID), 'FontSize', 14, 'FontWeight', 'bold');

    %% 9. Export Spreadsheets
    cd ../SensAnalytical/
    
    LSAoutputFolder = 'PatientDataLSA';
    
    if ~exist(LSAoutputFolder,'dir')
        mkdir(LSAoutputFolder);
    end

    comparisonFile = fullfile(LSAoutputFolder,'LSA_MethodsComparison.xlsx');
    stateFile      = fullfile(LSAoutputFolder,'LSA_AnalyticalStates.xlsx');
    
    param_labels = {'Rpul','Rmval','Raval','Rsys','Rtval','Rpval',...
                    'Cpv','Csa','Csv','Cpa','ElvM','Elvm','ErvM','Ervm','TC','TR'};
    
    % Build one row for this patient
    
    row = table(PID,'VariableNames',{'Patient'});
    
    for i = 1:num_params
    
        row.(sprintf('%s_Analytical',param_labels{i})) = analytical_overall(i);
    
        row.(sprintf('%s_FiniteDiff',param_labels{i})) = sens_norm_scaled(i);
    
        row.(sprintf('%s_AbsError',param_labels{i})) = abs_error(i);
    
        row.(sprintf('%s_RelErrorPct',param_labels{i})) = rel_error(i);
    
    end
    
    % Update existing spreadsheet if it exists
    % Export spreadsheet containing all patient sensitivities
    if isfile(comparisonFile)
    
        oldTable = readtable(comparisonFile);
    
        % Remove any existing row for this patient
        oldTable(oldTable.Patient==PID,:) = [];
    
        % Add updated row
        newTable = [oldTable; row];
    
    else
    
        newTable = row;
    
    end
    
    % Keep rows ordered by patient number
    newTable = sortrows(newTable,'Patient');
    
    % Save
    writetable(newTable,comparisonFile, 'UseExcel', false);
    
    fprintf('Updated %s\n',comparisonFile);
    
    % Export state-by-state analytical sensitivities
    
    state_names = {'Pulmonary Veins',...
                   'Left Ventricle',...
                   'Systemic Arteries',...
                   'Systemic Veins',...
                   'Right Ventricle',...
                   'Pulmonary Arteries'};
    
    stateRows = table();
    
    for s = 1:num_states
    
        temp = table(PID,string(state_names{s}),...
            'VariableNames',{'Patient','State'});
    
        for p = 1:num_params
            temp.(param_labels{p}) = Rsens_matrix(s,p);
        end
    
        stateRows = [stateRows; temp];
    
    end
    
    if isfile(stateFile)
    
        oldStates = readtable(stateFile);
    
        % Remove previous results for this patient
        oldStates(oldStates.Patient==PID,:) = [];
    
        stateTable = [oldStates; stateRows];
    
    else
    
        stateTable = stateRows;
    
    end
    
    stateTable = sortrows(stateTable,{'Patient','State'});
    
    writetable(stateTable,stateFile, 'UseExcel', false);
    
    fprintf('Updated %s\n',stateFile);

    %% 10. Automated Figure Saving Block
    outputFolder = 'PatientFiguresLSA';
    
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    
    imagePath1 = fullfile(outputFolder, sprintf('Patient_%d_CompartmentsRanked.png', PID));
    imagePath2 = fullfile(outputFolder, sprintf('Patient_%d_Ranked.png', PID));
    imagePath3 = fullfile(outputFolder, sprintf('Patient_%d_TimeSeries.png', PID));

    saveas(fig1, imagePath1);
    saveas(fig2, imagePath2);
    saveas(fig_time, imagePath3);
    
    fprintf('Figures successfully exported.');
    toc
end