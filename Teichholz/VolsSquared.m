function VolsSquared(PIDs)

close all;
% Comparing TTE Calculated Volume with MRI Volume Data

% SET UP AND CALCULATIONS
% Preallocate empty vectors
VolD_vec = []; % vector to store MRI max volumes
VolS_vec = []; % vector to store MRI min volumes
IDlvD_vec = []; % vector to store max inner diameter from TTE data
IDlvS_vec = []; % vector to store min inner diameter from TTE data

% badset = [82, 99, 152, 157, 198, 235, 281, 323]; % outliers

% change directory to extract patient data
cd ../Core/

for k = 1:length(PIDs)
    PID = PIDs(k);  
    % if ~(ismember(PID, badset))
        data = Patient(PID, 0);

        % Store MRI max and min volumes to vectors
        VolD_vec(k) = data.VlvM;
        VolS_vec(k) = data.VlvmT;
    
        % Store TTE max and min volumes to vectors
        IDlvD_vec(k) = data.IDlvD; % inner diameter at its largest during diastole, measured in cm
        IDlvS_vec(k) = data.IDlvS; % inner diameter at its largest during systole, measured in cm
    % end
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

% CALCULATING BASED ON V = D^2
% Parameter estimation
guess2 = [1,0];
pars2 = fminsearch(@(p) squared_error(p, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec), guess2);

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

figure(1); % V = D^2 scatter plot
% subplot(1,2,1);
    C2 = [ones(size(Vdias2)), 2*ones(size(Vsys2))];
    scatter(Vols2, Vols, 36, C2, 'filled');
    colormap([0 0 0.5; ...
            0.5 0.5 1]);
    hold on;
    y_fit2_quad = polyval(coeffs2_quad, Vols2_sorted);
    y_fit2_lin = polyval(coeffs2_lin, Vols2_sorted);
    reg2_quad = plot(Vols2_sorted, y_fit2_quad, 'b', 'LineWidth', 2, 'LineStyle', '-');
    reg2_lin = plot(Vols2_sorted, y_fit2_lin, 'r', 'LineWidth', 2, 'LineStyle', '-');
    d2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
    h2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold off;
    title('V = aD^2+b', 'FontSize', 15);
    xlabel('TTE Volume (mL)', 'FontSize', 15);
    ylabel('MRI Volume (mL)', 'FontSize', 15);
    legend([d2, h2, reg2_quad, reg2_lin], ...
           {'LV Dias', 'LV Sys', 'Quadratic Regression', 'Linear Regression'}, ...
           'Location', 'southeast', 'FontSize', 15);
    regText = sprintf('R^2 for Linear = %.4f\nR^2 for Quadratic = %.4f', R2_lin2, R2_quad2);
    slopeText = sprintf('y = x\ny = %.4fx^2+%.4fx+%.4f', coeffs2_quad(1), coeffs2_quad(2), coeffs2_quad(3));
    annotation('textbox', [.15 .7 .15 .13], 'String', regText, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 15);
    annotation('textbox', [.15 .6 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 15);
figure(2); % V = D^2
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

disp('Displaying Slopes for V = aD^2 + b')
coeffs2_lin
coeffs2_quad
coeffs2_log

% CONFIDENCE INTERVAL AND TEST STATISTICS
[h, p, ci, stats] = ttest2(Vols, Vols2);
disp('95% Confidence Interval for the difference:');
disp(ci);

error = rmse(Vols2, Vols) %rmse

% Proportional bias check
p = polyfit(avg2, diff2, 1);
slope = p(1); 
intercept = p(2);

[r, pval] = corr(avg2', diff2');

disp('Proportional Bias Check')
fprintf('Slope: %.4f\n', slope)
fprintf('Intercept: %.4f\n', intercept)
fprintf('Correlation (r): %.4f\n', r)
fprintf('p-value: %.4f\n', pval)

if pval < 0.05
    disp('Significant proportioanl bias detected. Differences change systematically with volume magnitude.')
else
    disp('No significant proportioanl bias. Differences appear consistent across the volume range.')
end

figure(3);
scatter(avg2, diff2, 'filled');
hold on;
yline(bias2, 'k-', 'LineWidth', 1.5);
yline(uLOA2, 'r--', 'LineWidth', 1.2);
yline(lLOA2, 'r--', 'LineWidth', 1.2);
xFit = linspace(min(avg2), max(avg2), 100);
yFit = polyval(p, xFit);
plot(xFit,yFit, 'b-', 'LineWidth', 1.5);
xlabel('Mean of TTE and MRI Volumes')
ylabel('Difference (MRI - TTE)');
title('Bland-Altman Plot with Proportioanl Bias Line');
legend('Data points', 'Mean difference', 'Upper LOA', 'Lower LOA', 'Bias trend line', 'Location', 'best');
hold off;

end

%%  Nested error function for V = D^2
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