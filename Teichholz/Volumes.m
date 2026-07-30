function Volumes(PIDs)
close all;
% Comparing TTE Calculated Volume with MRI Volume Data

% SET UP AND CALCULATIONS
% Preallocate empty vectors
VolD_vec = []; % vector to store MRI max volumes
VolS_vec = []; % vector to store MRI min volumes
IDlvD_vec = []; % vector to store max inner diameter from TTE data
IDlvS_vec = []; % vector to store min inner diameter from TTE data

badset = [82, 99, 152, 157, 198, 235, 281, 323]; % outliers

% change directory to extract patient data
cd ../Core/

for k = 1:length(PIDs)
    PID = PIDs(k);  
    if ~(ismember(PID, badset))
        data = Patient(PID, 0);

        % Store MRI max and min volumes to vectors
        VolD_vec(k) = data.VlvM;
        VolS_vec(k) = data.VlvmT;
    
        % Store TTE max and min volumes to vectors
        IDlvD_vec(k) = data.IDlvD; % inner diameter at its largest during diastole, measured in cm
        IDlvS_vec(k) = data.IDlvS; % inner diameter at its largest during systole, measured in cm
    end
end

% back to this directory
cd ../Teichholz/

% Check that vector entries are okay
good = ~isnan(VolD_vec) & ~isnan(VolS_vec) & ~isnan(IDlvD_vec) & ~isnan(IDlvS_vec) ...
       & VolD_vec > 0 & VolS_vec > 0 & IDlvD_vec > 0 & IDlvS_vec > 0;
VolD_vec = VolD_vec(good);
VolS_vec = VolS_vec(good);
IDlvD_vec = IDlvD_vec(good);
IDlvS_vec = IDlvS_vec(good);

% Combine vectors
Vols = [VolD_vec, VolS_vec];
IDlv = [IDlvD_vec, IDlvS_vec];

% Calculation of volume
VolT = (7.*IDlv.^3)./(2.4+IDlv);
% VolT = (pi.*IDlv.^3)./(6.*((0.075.*IDlv)+ 0.18));

% INFORMATION FOR PLOTS
% Linear regression
coeffsT = polyfit(Vols, VolT, 1);
mT = coeffsT(1);

% To plot y = x
xRange = linspace(min(Vols), max(Vols), 100);

% Bland-Altman Plot Calculations
avgT = (Vols+VolT)/2;
diffT = Vols - VolT;
biasT = mean(diffT);
sdDiffT = std(diffT);
uLOAT = biasT + (1.96*sdDiffT);
lLOAT = biasT - (1.96*sdDiffT);

%% CALCULATING BASED ON V = D^3
% Parameter estimation
guess = [1.2, 0];
pars = fminsearch(@(p) errors(p, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec), guess);

% sse = errors(pars, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec);
% disp(sse)
% pause;

pars_a = pars(1);
pars_b = pars(2);
disp('Parameters for V = aD^3 + b')
disp(pars_a)
disp(pars_b)

% Calculations for volume
Vdias = pars_a.*IDlvD_vec.^3 + pars_b;
Vsys = pars_a.*IDlvS_vec.^3 + pars_b;

% Regression
Vols3 = [Vdias, Vsys];
coeffs3_log = polyfit(log(Vols3), Vols, 1); % logarithmic
coeffs3_quad = polyfit(Vols3, Vols, 2); % quadratic
coeffs3_lin = polyfit(Vols3, Vols, 1);
% coeffs3_power = polyfit(log(Vols3), log(Vols), 1); % power
% coeffs3_exp = polyfit(Vols3, log(Vols), 1); % exponential

% Sort to plot
[Vols3_sorted, sortIdx] = sort(Vols3);

% Bland-Altman Plot information
avg3 = (Vols+Vols3)/2;
diff3 = Vols - Vols3;
bias3 = mean(diff3);
sdDiff3 = std(diff3);
uLOA3 = bias3 + (1.96*sdDiff3);
lLOA3 = bias3 - (1.96*sdDiff3);

% R^2 Calculations
disp('R^2 for V = aD^3 + b')
mdl_linear3 = fitlm(Vols3, Vols);
R2_lin3 = mdl_linear3.Rsquared.Ordinary;
disp('Coefficient of Determination for Linear Fit:')
disp(R2_lin3)

mdl_quad3 = fitlm(Vols3, Vols, 'purequadratic');
R2_quad3 = mdl_quad3.Rsquared.Ordinary;
disp('Coefficient of Determination for Quadratic Fit:')
disp(R2_quad3)

mdl_log3 = fitlm(log(Vols3), Vols);
R2_log3 = mdl_log3.Rsquared.Ordinary;
disp('Coefficient of Determination for Logarithmic Fit:')
disp(R2_log3)

% mdl_exp3 = fitlm(Vols3, log(Vols));
% R2_exp3 = mdl_exp3.Rsquared.Ordinary;
% disp('Coefficient of Determination for Exponential Fit:')
% disp(R2_exp3)
% 
% mdl_power3 = fitlm(log(Vols3), log(Vols));
% R2_power3 = mdl_power3.Rsquared.Ordinary;
% disp('Coefficient of Determination for Power Fit')
% disp(R2_power3)



%% CALCULATING BASED ON V = D^2
% Parameter estimation
guess2 = [1,0];
pars2 = fminsearch(@(p) squared_error(p, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec), guess2);

% sse2 = squared_error(pars2, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec);

pars2_a = pars2(1);
pars2_b = pars2(2);
disp('Parameters for V = aD^2 + b')
disp(pars2_a)
disp(pars2_b)

% Calculations for volume
Vdias2 = pars2_a.*IDlvD_vec.^2 + pars2_b;
Vsys2 = pars2_a.*IDlvS_vec.^2 + pars2_b;

% Regression
Vols2 = [Vdias2, Vsys2];
coeffs2_log = polyfit(log(Vols2), Vols, 1);
coeffs2_lin = polyfit(Vols2, Vols, 1);
coeffs2_quad = polyfit(Vols2, Vols, 2);
% coeffs2_power = polyfit(log(Vols2), log(Vols), 1); % power
% coeffs2_exp = polyfit(Vols2, log(Vols), 1); % exponential

% Sort to plot
[Vols2_sorted, sortIdx2] = sort(Vols2);

% Bland-Altman Plot information
avg2 = (Vols+Vols2)/2;
diff2 = Vols - Vols2;
bias2 = mean(diff2);
sdDiff2 = std(diff2);
uLOA2 = bias2 + (1.96*sdDiff2);
lLOA2 = bias2 - (1.96*sdDiff2);

% R^2 Calculations
disp('R^2 for V = aD^2 + b')
mdl_linear2 = fitlm(Vols2, Vols);
R2_lin2 = mdl_linear2.Rsquared.Ordinary;
disp('Coefficient of Determination for Linear Fit:')
disp(R2_lin2)

mdl_quad2 = fitlm(Vols2, Vols, 'purequadratic');
R2_quad2 = mdl_quad2.Rsquared.Ordinary;
disp('Coefficient of Determination for Quadratic Fit:')
disp(R2_quad2)

mdl_log2 = fitlm(log(Vols2), Vols);
R2_log2 = mdl_log2.Rsquared.Ordinary;
disp('Coefficient of Determination for Logathrimic Fit:')
disp(R2_log2)

% mdl_exp2 = fitlm(Vols2, log(Vols));
% R2_exp2 = mdl_exp2.Rsquared.Ordinary;
% disp('Coefficient of Determination for Exponential Fit:')
% disp(R2_exp2)
% 
% mdl_power2 = fitlm(log(Vols2), log(Vols));
% R2_power2 = mdl_power2.Rsquared.Ordinary;
% disp('Coefficient of Determination for Power Fit')
% disp(R2_power2)

%% CALCULATING BASED ON V = a*pi*r^3
guess = 1;
parspi = fminsearch(@(p) pi_error(p, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec), guess);

pi_a = parspi(1);
disp('Parameter for V = a*pi*r^3:')
disp(pi_a)

% Calculations for volume
R_dias = IDlvD_vec./2;
R_sys = IDlvS_vec./2;
Vdias_pi = pi_a.*pi.*R_dias.^3;
Vsys_pi = pi_a.*pi.*R_sys.^3;

% Regression
Volspi = [Vdias_pi, Vsys_pi];

coeffspi_lin = polyfit(Volspi, Vols, 1);
coeffspi_log = polyfit(log(Volspi), Vols, 1);
coeffspi_quad = polyfit(Volspi, Vols, 2);
% coeffspi_power = polyfit(log(Volspi), log(Vols), 1);
% coeffspi_exp = polyfit(Volspi, log(Vols), 1);

% sorting to plot
[Volspi_sorted, sortIdxpi] = sort(Volspi);

% Bland-Altman Plot information
avgpi = (Volspi+Vols)/2;
diffpi = Vols - Volspi;
biaspi = mean(diffpi);
sdDiffpi = std(diffpi);
uLOApi = biaspi + (1.96*sdDiffpi);
lLOApi = biaspi - (1.96*sdDiffpi);

% R^2 Calculations
disp('R^2 for V = a*pi*r^3')
mdl_lin_pi = fitlm(Volspi, Vols);
R2_linpi = mdl_lin_pi.Rsquared.Ordinary;
disp('Coefficient of Determination for Linear Fit:')
disp(R2_linpi)

mdl_quadpi = fitlm(Volspi, Vols, 'purequadratic');
R2_quadpi = mdl_quadpi.Rsquared.Ordinary;
disp('Coefficient of Determination for Quadratic Fit:')
disp(R2_quadpi)

mdl_logpi = fitlm(log(Volspi), Vols);
R2_logpi = mdl_logpi.Rsquared.Ordinary;
disp('Coefficient of Determination for Logathrimic Fit:')
disp(R2_logpi)

% mdl_exppi = fitlm(Volspi, log(Vols));
% R2_exppi = mdl_exppi.Rsquared.Ordinary;
% disp('Coefficient of Determination for Exponential Fit:')
% disp(R2_exppi)

% mdl_powerpi = fitlm(log(Volspi), log(Vols));
% R2_powerpi = mdl_powerpi.Rsquared.Ordinary;
% disp('Coefficient of Determination for Power Fit')
% disp(R2_powerpi)

%% PLOTS
figure(1); % Teichholz equation
subplot(1,3,1); % scatter plot
    scatter(Vols, VolT, 'filled', 'MarkerFaceColor','b'); % data and calculated values
    hold on;
    y_fitT = polyval(coeffsT, Vols);
    plot(Vols, y_fitT, 'r', 'LineWidth', 2, 'LineStyle', '-'); % regression
    plot(xRange, xRange, 'm', 'LineWidth',2, 'LineStyle', ':'); % y = x
    hold off;
    title('MRI vs TTE Left Ventricular Volume at Diastole');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend('Data Points', 'Regression Line', 'y=x','Location', 'southeast');
    slopeText = sprintf('Slope = %.2f', mT);
    annotation('textbox', [.15 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
subplot(1,3,2); % bland-altman plot
    plot(avgT, diffT, 'ko', 'MarkerFaceColor','b', 'MarkerSize',8);
    hold on;
    xlims = xlim;
    line(xlims, [biasT, biasT], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOAT, uLOAT], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOAT, lLOAT], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasT + (max(diffT)*0.05), sprintf('Bias: %.2f', biasT), 'Color', 'r');
    text(xlims(2)*0.75, uLOAT + (max(diffT)*0.05), sprintf('+1.96 SD: %.2f', uLOAT), 'Color', 'r');
    text(xlims(2)*0.75, lLOAT - (max(diffT)*0.05), sprintf('-1.96 SD: %.2f', lLOAT), 'Color', 'r');
    grid on;
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE (MRI - TTE)');
    title('Bland-Altman Plot for Diastolic Volume of LV');
    hold off;
subplot(1,3,3); % box and whiskers plot
    bT = boxplot(diffT);
    ylabel('Difference Between MRI and TTE (MRI - TTE)');
    set(bT, 'linewidth', 2);
figure(2); % V = D^3
C3 = [ones(size(Vdias)), 2*ones(size(Vsys))];
% scatter plot
    scatter(Vols3, Vols, 36, C3, 'filled');
    colormap([0 0 0.5; ...
            0.5 0.5 1]);
    hold on;
    % y_fit3 = polyval(coeffs3, log(Vols3));
    % y_fit3_log = polyval(coeffs3_log, log(Vols3_sorted));
    y_fit3_quad = polyval(coeffs3_quad, Vols3_sorted);
    y_fit3_lin = polyval(coeffs3_lin, Vols3_sorted);
    % y_fit3_power = polyval(coeffs3_power, log(Vols3_sorted));
    % y_fit3_exp = polyval(coeffs3_exp, Vols3_sorted);
    % reg3_log = plot(Vols3_sorted, y_fit3_log, 'r', 'LineWidth',2, 'LineStyle', '-');
    reg3_quad = plot(Vols3_sorted, y_fit3_quad, 'b', 'LineWidth', 2, 'LineStyle', '-');
    reg3_lin = plot(Vols3_sorted, y_fit3_lin, 'Color', [0 0.5 0], 'LineWidth', 2, 'LineStyle', '-');
    % reg3_power = plot(Vols3_sorted, y_fit3_power, 'Color', [0.92, 0.65, 0.0], 'LineWidth', 2, 'LineStyle', '-');
    % reg3_exp = plot(Vols3_sorted, y_fit3_exp, 'Color', [0.65, 0.0, 0.45], 'LineWidth', 2, 'LineStyle', '-');
    d3 = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
    h3 = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold off;
    title('Calculating Volume as V = aD^3+b');
    xlabel('TTE Volume (mL)');
    ylabel('MRI Volume (mL)');
    legend([d3, h3, reg3_quad, reg3_lin], ...
           {'LV Dias', 'LV Sys','Quadratic Regression', 'Linear Regression'}, ...
           'Location', 'southeast');  
    regText = sprintf('R^2 for Linear = %.4f\nR^2 for Quadratic = %.4f', R2_lin3, R2_quad3);
    annotation('textbox', [.15 .7 .15 .13], 'String', regText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    % slopeText = sprintf('Slope = %.2f', m3);
    % annotation('textbox', [.15 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    % slopeText = sprintf('Linear Slope = %.4\nQuadratic Slope = %.4', coeffs3_lin(1), coeffs3_quad(1));
    % annotation('textbox', [.25 .17 .25 .23], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
figure(3); 
subplot(1,2,1)% Bland-Altman plot
    scatter(avg3, diff3, 50, C3, 'filled', 'MarkerEdgeColor', 'k');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    xlims = xlim;
    line(xlims, [bias3, bias3], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOA3, uLOA3], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOA3, lLOA3], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, bias3 + (max(diff3)*0.05), sprintf('Bias: %.2f', bias3), 'Color', 'r');
    text(xlims(2)*0.75, uLOA3 + (max(diff3)*0.05), sprintf('+1.96 SD: %.2f', uLOA3), 'Color', 'r');
    text(xlims(2)*0.75, lLOA3 - (max(diff3)*0.05), sprintf('-1.96 SD: %.2f', lLOA3), 'Color', 'r');
    grid on;
    pd = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0 0 0.5]);
    ps = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0.5 0.5 1]);
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE');
    title('Bland-Altman Plot for V = aD^3+b');
    legend([pd, ps], ...
           {'LV Dias', 'LV Sys'}, ...
           'Location','southwest');
    hold off;
subplot(1,2,2); % box and whiskers plot
    b3 = boxplot(diff3);
    title('Box and Whiskers Plot for V = aD^3+b');
    set(b3, 'LineWidth', 2);

figure(4); % V = D^2 scatter plot
    C2 = [ones(size(Vdias2)), 2*ones(size(Vsys2))];
    scatter(Vols2, Vols, 36, C2, 'filled');
    colormap([0 0 0.5; ...
            0.5 0.5 1]);
    hold on;
    % y_fit3 = polyval(coeffs3, log(Vols3));
    % y_fit2_log = polyval(coeffs2_log, log(Vols2_sorted));
    y_fit2_quad = polyval(coeffs2_quad, Vols2_sorted);
    y_fit2_lin = polyval(coeffs2_lin, Vols2_sorted);
    % y_fit2_power = polyval(coeffs2_power, log(Vols2_sorted));
    % y_fit2_exp = polyval(coeffs2_exp, Vols2_sorted);
    % reg2_log = plot(Vols2_sorted, y_fit2_log, 'r', 'LineWidth',2, 'LineStyle', '-');
    reg2_quad = plot(Vols2_sorted, y_fit2_quad, 'b', 'LineWidth', 2, 'LineStyle', '-');
    reg2_lin = plot(Vols2_sorted, y_fit2_lin, 'Color', [0 0.5 0], 'LineWidth', 2, 'LineStyle', '-');
    % reg2_power = plot(Vols2_sorted, y_fit2_power, 'Color', [0.92, 0.65, 0.0], 'LineWidth', 2, 'LineStyle', '-');
    % reg2_exp = plot(Vols2_sorted, y_fit2_exp, 'Color', [0.65, 0.0, 0.45], 'LineWidth', 2, 'LineStyle', '-');
    d2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
    h2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold off;
    title('Calculating Volume as V = aD^2+b');
    xlabel('TTE Volume (mL)');
    ylabel('MRI Volume (mL)');
    legend([d2, h2, reg2_quad, reg2_lin], ...
           {'LV Dias', 'LV Sys', 'Quadratic Regression', 'Linear Regression'}, ...
           'Location', 'southeast');
    regText = sprintf('R^2 for Linear = %.4f\nR^2 for Quadratic = %.4f', R2_lin2, R2_quad2);
    annotation('textbox', [.15 .7 .15 .13], 'String', regText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    % slopeText = sprintf('Slope = %.2f', m3);
    % annotation('textbox', [.15 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    % slopeText = sprintf('Linear Slope = %.4\nQuadratic Slope = %.4', coeffs2_lin(1), coeffs2_quad(1));
    % annotation('textbox', [.25 .17 .25 .23], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
figure(5); % V = D^2
subplot(1,2,1); % bland-altman plot
    scatter(avg2, diff2, 50, C2, 'filled', 'MarkerEdgeColor', 'k');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    xlims = xlim;
    line(xlims, [bias2, bias2], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOA2, uLOA2], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOA2, lLOA2], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, bias2 + (max(diff2)*0.05), sprintf('Bias: %.2f', bias2), 'Color', 'r');
    text(xlims(2)*0.75, uLOA2 + (max(diff2)*0.05), sprintf('+1.96 SD: %.2f', uLOA2), 'Color', 'r');
    text(xlims(2)*0.75, lLOA2 - (max(diff2)*0.05), sprintf('-1.96 SD: %.2f', lLOA2), 'Color', 'r');
    grid on;
    pd = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0 0 0.5]);
    ps = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0.5 0.5 1]);
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE');
    title('Bland-Altman Plot for V = aD^2+b');
    legend([pd, ps], ...
           {'LV Dias', 'LV Sys'}, ...
           'Location','southwest');
    hold off;
subplot(1,2,2); % box and whisker plot
    b2 = boxplot(diff2);
    title('Box and Whiskers Plot for V = aD^2+b');
    set(b2, 'LineWidth', 2);

figure(6); % scatter plot for V = a*pi*R^3
Cpi = [ones(size(Vdias_pi)), 2*ones(size(Vsys_pi))];
    scatter(Volspi, Vols, 36, Cpi, 'filled');
    colormap([0 0 0.5; ...
            0.5 0.5 1]);
    hold on;
    y_fitpi_log = polyval(coeffspi_log, log(Volspi_sorted));
    y_fitpi_quad = polyval(coeffspi_quad, Volspi_sorted);
    y_fitpi_lin = polyval(coeffspi_lin, Volspi_sorted);
    % y_fitpi_power = polyval(coeffspi_power, log(Volspi_sorted));
    % y_fitpi_exp = polyval(coeffspi_exp, Volspi_sorted);
    regpi_log = plot(Volspi_sorted, y_fitpi_log, 'r', 'LineWidth',2, 'LineStyle', '-');
    regpi_quad = plot(Volspi_sorted, y_fitpi_quad, 'b', 'LineWidth', 2, 'LineStyle', '-');
    regpi_lin = plot(Volspi_sorted, y_fitpi_lin, 'Color', [0 0.5 0], 'LineWidth', 2, 'LineStyle', '-');
    % regpi_power = plot(Volspi_sorted, y_fitpi_power, 'Color', [0.92, 0.65, 0.0], 'LineWidth', 2, 'LineStyle', '-');
    % regpi_exp = plot(Volspi_sorted, y_fitpi_exp, 'Color', [0.65, 0.0, 0.45], 'LineWidth', 2, 'LineStyle', '-');
    dpi = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
    hpi = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold off;
    title('Calculating Volume as V = a\piR^3');
    xlabel('TTE Volume (mL)');
    ylabel('MRI Volume (mL)');
    legend([dpi, hpi, regpi_lin, regpi_quad, regpi_log], ...
           {'LV Dias', 'LV Sys', 'Linear Regression', 'Quadratic Regression', 'Logarithmic Regression'}, ...
           'Location', 'southeast');
    regText = sprintf('R^2 for Linear = %.4f\nR^2 for Quadratic = %.4f\nR^2 for Logarithmic = %.4f', ...
        R2_linpi, R2_quadpi, R2_logpi);
    annotation('textbox', [.15 .7 .15 .13], 'String', regText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    % slopeText = sprintf('Linear Slope = %.4\nLogarithmic Slope = %.4\nQuadratic Slope = %.4', coeffspi_lin(1), coeffspi_log(1), coeffspi_quad(1));
    % annotation('textbox', [.25 .17 .25 .23], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
figure(7); % V = a*pi*R^3
subplot(1,2,1); % bland-altman plot
    scatter(avgpi, diffpi, 50, Cpi, 'filled', 'MarkerEdgeColor', 'k');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    xlims = xlim;
    line(xlims, [biaspi, biaspi], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOApi, uLOApi], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOApi, lLOApi], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biaspi + (max(diffpi)*0.05), sprintf('Bias: %.2f', biaspi), 'Color', 'r');
    text(xlims(2)*0.75, uLOApi + (max(diffpi)*0.05), sprintf('+1.96 SD: %.2f', uLOApi), 'Color', 'r');
    text(xlims(2)*0.75, lLOApi - (max(diffpi)*0.05), sprintf('-1.96 SD: %.2f', lLOApi), 'Color', 'r');
    grid on;
    pd = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0 0 0.5]);
    ps = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0.5 0.5 1]);
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE');
    title('Bland-Altman Plot for V = a\piR^3');
    legend([pd, ps], ...
           {'LV Dias', 'LV Sys'}, ...
           'Location','southwest');
    hold off;
subplot(1,2,2); % box and whisker plot
    bpi = boxplot(diffpi);
    title('Box and Whiskers Plot for V = a\piR^3');
    set(bpi, 'LineWidth', 2);

%% Displaying Slopes
disp('Displaying Slopes for V = aD^3 + b')
coeffs3_lin
coeffs3_quad
coeffs3_log

disp('Displaying Slopes for V = aD^2 + b')
coeffs2_lin
coeffs2_quad
coeffs2_log

disp('Displaying Slopes for V = apiR^3')
coeffspi_lin
coeffspi_quad
coeffspi_log


%% CONFIDENCE INTERVAL
[h, p, ci, stats] = ttest2(Vols, Vols2, 'Alpha', 0.05);
disp('95% Confidence Interval for the difference:');
disp(ci);

end

%% Nested Error Functions
% Nested error function for V = D^3
function sse = errors(params, IDlvD, IDlvS, VolD_real, VolS_real)
    a = params(1);
    b = params(2);

    % Combine vectors
    Vol_real = [VolD_real, VolS_real];
    IDlv = [IDlvD, IDlvS];

    % Volume calculation
    Vol_calc = a.*IDlv.^3 + b;

    % Sum of squared errors to minimize distance between calculations
    sse = sum((Vol_real - Vol_calc).^2);
end

% Nested error function for V = D^2
function sse2 = squared_error(params, IDlvD, IDlvS, VolD_real, VolS_real)
    a = params(1);
    b = params(2);

    % Combine vectors
    Vol_real = [VolD_real, VolS_real];
    IDlv = [IDlvD, IDlvS];

    % Volume calculation
    Vol_calc = a.*IDlv.^2 + b;

    % Sum of squared errors to minimize distance between calculations
    sse2 = sum((Vol_real - Vol_calc).^2);
end

% Nested error function for V = a*pi*r^3
function ssepi = pi_error(params, IDlvD, IDlvS, VolD_real, VolS_real)
    a = params(1);

    % Calculate radius
    RD = IDlvD./2;
    RS = IDlvS./2;

    % Combine vectors
    Vold_real = [VolD_real, VolS_real];
    R = [RD, RS];

    % Volume Calculations
    Vol_calc = a.*pi.*R.^3;

    % Sum of squared errors to minimize distance between calculations
    ssepi = sum((Vold_real - Vol_calc).^2);
end