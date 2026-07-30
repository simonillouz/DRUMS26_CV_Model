function Scattered_Vols(PIDs)
close all;

% PREALLOCATED EMPTY VECTORS
VolD_vec = []; % vector to store MRI max volumes
VolS_vec = []; % vector to store MRI min volumes
IDlvD_vec = []; % vector to store TTE diastolic volumes
IDlvS_vec = []; % vector to store TTE systolic volumes

% change directory to extract patient data
cd ../Core/
for k = 1:length(PIDs)
    PID = PIDs(k);
    data = Patient(PID);

    % Store MRI max and min volumes to vectors
    VolD_vec(k) = data.VlvM;
    VolS_vec(k) = data.VlvmT; 

    % Store TTE max and min volumes to vectors
    IDlvD_vec(k) = data.IDlvD;
    IDlvS_vec(k) = data.IDlvS;
end

% back to this directory
cd ../Teichholz/

% Check that vectors are okay
good = ~isnan(VolD_vec) & ~isnan(VolS_vec) & ~isnan(IDlvD_vec) & ~isnan(IDlvS_vec);
VolD_vec = VolD_vec(good);
VolS_vec = VolS_vec(good);
IDlvD_vec = IDlvD_vec(good);
IDlvS_vec = IDlvS_vec(good);

% Initial guess for a and b (used for both diastolic and systoliv volume)
guess = [0.0425, 0.125];

% PARAMETER OPTIMIZATION FOR A AND B - MINIMIZING DISTANCE FROM SLOPE OF 1
parms_slope = fminsearch(@(p) error_slope(p, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec), guess);
a_slope = parms_slope(1);
b_slope = parms_slope(2);
disp('Parameters for Minimizig Distance from Slope:')
disp(a_slope)
disp(b_slope)

% MINIMIZING DISTANCE FROM ACTUAL VALUES
parms_val = fminsearch(@(p) error_val(p, IDlvD_vec, IDlvS_vec, VolD_vec, VolS_vec), guess);
a_val = parms_val(1);
b_val = parms_val(2);
disp('Parameters for Minimizing Distance from Values:')
disp(a_val)
disp(b_val)

% CALCULATIONS FOR VOLUME - SLOPE
% Diastolic
LD = 1./(a_slope*IDlvD_vec+b_slope);
VolDT_slope = 0.5.*((pi/3).*(IDlvD_vec.^3).*LD);
% Systolic
LS = 1./(a_slope*IDlvS_vec+b_slope);
VolST_slope = 0.5.*((pi/3).*(IDlvS_vec.^3).*LS);

% CALCULATIONS FOR VOLUME - VALUES
% Diastolic
DD = 1./(a_val*IDlvD_vec+b_val);
VolDT_val = 0.5.*((pi/3).*(IDlvD_vec.^3).*DD);
% Systolic
DS = 1./(a_val*IDlvS_vec+b_val);
VolST_val = 0.5.*((pi/3).*(IDlvS_vec.^3).*DS);

% LINEAR REGRESSION - SLOPE
Vols = [VolD_vec, VolS_vec];
VolT_slope = [VolDT_slope, VolST_slope];
coeffs_slope = polyfit(Vols, VolT_slope, 1);
m_slope = coeffs_slope(1);

% Calculations of R^2 - Slope
mdl_slope_lin = fitlm(Vols, VolT_slope);
R2_linslope = mdl_slope_lin.Rsquared.Ordinary;
disp('Slope - Coefficient of Determination for Linear Fit:')
disp(R2_linslope)

% LINEAR REGRESSION - VALUES
VolT_val = [VolDT_val, VolST_val];
coeffs_val = polyfit(Vols, VolT_val, 1);
m_val = coeffs_val(1);

% Calculations of R^2 - Values
mdl_val_lin = fitlm(Vols, VolT_val);
R2_linval = mdl_val_lin.Rsquared.Ordinary;
disp('Values - Coefficient of Determination of Linear Fit:')
disp(R2_linval)

% FOR GRAPHS
xRange = linspace(min(Vols), max(Vols), 100);

avg_slope = (Vols+VolT_slope)/2;
avg_val = (Vols+VolT_val)/2;

diff_slope = Vols - VolT_slope;
diff_val = Vols - VolT_val;

bias_slope = mean(diff_slope);
bias_val = mean(diff_val);

sdDiff_slope = std(diff_slope);
sdDiff_val = std(diff_val);

uLOA_slope = bias_slope + (1.96*sdDiff_slope);
lLOA_slope = bias_slope - (1.96*sdDiff_slope);
uLOA_val = bias_val + (1.96*sdDiff_val);
lLOA_val = bias_val - (1.96*sdDiff_val);

% PLOTS
% To differentiate colors
C_slope = [ones(size(VolDT_slope)), 2*ones(size(VolST_slope))];
C_val = [ones(size(VolDT_val)), 2*ones(size(VolST_val))];

figure(1); % scatter plots
subplot(1,2,1); % minimizing slope
    scatter(Vols, VolT_slope, 36, C_slope, 'filled');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    y_fit_slope = polyval(coeffs_slope, Vols);
    hreg = plot(Vols, y_fit_slope, 'r','LineWidth', 2, 'LineStyle', '-'); % regression
    hid = plot(xRange, xRange, 'm', 'LineWidth',2, 'LineStyle', ':'); % y = x
    hd = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
    hs = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold off;
    title('Optimization with Minimizing Slope');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend([hd, hs, hreg, hid], ...
           {'LV Dias', 'LV Sys', 'Regression Line', 'y=x'}, ...
           'Location', 'southeast');    
    slopeText = sprintf('Slope = %.2f', m_slope);
    annotation('textbox', [.15 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
subplot(1,2,2); % minimizing measurement
    scatter(Vols, VolT_val, 36, C_val, 'filled');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    y_fit_val = polyval(coeffs_val, Vols);
    greg = plot(Vols, y_fit_val, 'r','LineWidth', 2, 'LineStyle', '-'); % regression
    gid = plot(xRange, xRange, 'm', 'LineWidth',2, 'LineStyle', ':'); % y = x
    gd = plot(NaN, NaN, 'o', 'MarkerFaceColor',[0 0 0.5], 'MarkerEdgeColor', 'none');
    gs = plot(NaN, NaN, 'o', 'MarkerFaceColor',[0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold off;
    title('Optimization with Minimizing Measurements');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend([gd, gs, greg, gid], ...
           {'LV Dias', 'LV Sys', 'Regression Line', 'y=x'}, ...
           'Location', 'southeast');       
    slopeText = sprintf('Slope = %.2f', m_val);
    annotation('textbox', [.60 .7 .1 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
figure(2); % bland-altman plots
subplot(1,2,1); % minimizing slope
    scatter(avg_slope, diff_slope, 50, C_slope, 'filled', 'MarkerEdgeColor', 'k');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    xlims = xlim;
    line(xlims, [bias_slope, bias_slope], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOA_slope, uLOA_slope], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOA_slope, lLOA_slope], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, bias_slope + (max(diff_slope)*0.05), sprintf('Bias: %.2f', bias_slope), 'Color', 'r');
    text(xlims(2)*0.75, uLOA_slope + (max(diff_slope)*0.05), sprintf('+1.96 SD: %.2f', uLOA_slope), 'Color', 'r');
    text(xlims(2)*0.75, lLOA_slope - (max(diff_slope)*0.05), sprintf('-1.96 SD: %.2f', lLOA_slope), 'Color', 'r');
    grid on;
    pd = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0 0 0.5]);
    ps = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0.5 0.5 1]);
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE');
    title('Minimizing Slope');
    legend([pd, ps], ...
           {'LV Dias', 'LV Sys'}, ...
           'Location','southeast');
    hold off;
subplot(1,2,2); % minimizing measurement
    scatter(avg_val, diff_val, 50, C_val, 'filled', 'MarkerEdgeColor', 'k');
    colormap([0 0 0.5; ... 
          0.5 0.5 1]);
    hold on;
    xlims = xlim;
    line(xlims, [bias_val, bias_val], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOA_val, uLOA_val], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOA_val, lLOA_val], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, bias_val + (max(diff_val)*0.05), sprintf('Bias: %.2f', bias_val), 'Color', 'r');
    text(xlims(2)*0.75, uLOA_val + (max(diff_val)*0.05), sprintf('+1.96 SD: %.2f', uLOA_val), 'Color', 'r');
    text(xlims(2)*0.75, lLOA_val - (max(diff_val)*0.05), sprintf('-1.96 SD: %.2f', lLOA_val), 'Color', 'r');
    grid on;
    md = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0 0 0.5]);
    cs = plot(NaN, NaN, 'ko', 'MarkerFaceColor', [0.5 0.5 1]);
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE');
    title('Minimizing Measurement');
    legend([md, cs], ...
           {'LV Dias', 'LV Sys'}, ...
           'Location','southeast');
    hold off;
figure(3); % box and whisker plots
subplot(1,2,1); % minimizing slope
    b = boxplot(diff_slope);
    ylabel('Difference Between MRI and TTE');
    title('Minimzing Slope');
    set(b, 'LineWidth', 2);
subplot(1,2,2); % minimizing measurement
    h = boxplot(diff_val);
    title('Minimzing Measurements');
    set(h, 'LineWidth', 2);
end
%% NESTED SUBFUNCTIONS
% Error function for minimizing distance of slope from 1
function sse_slope = error_slope(params, IDlvD, IDlvS, VolD_real, VolS_real)
    a = params(1);
    b = params(2);

    % Combine vectors
    Vol_real = [VolD_real, VolS_real];
    IDlv = [IDlvD, IDlvS];

    % Volume calculation
    L = 1./(a.*IDlv+b);
    VolT_calc = 0.5.*((pi/3).*(IDlv.^3).*L);

    % Calculate regression 
    coeffs = polyfit(Vol_real, VolT_calc, 1);
    m = coeffs(1);

    % Sum of squared errors to minimize distance between slopes
    sse_slope = (m-1)^2;
end

% Error function for minimizing distance between actual and calculated values
function sse_val = error_val(params, IDlvD, IDlvS, VolD_real, VolS_real)
    a = params(1);
    b = params(2);

    % Combine vectors
    Vol_real = [VolD_real, VolS_real];
    IDlv = [IDlvD, IDlvS];

    % Volume calculation
    L = 1./(a.*IDlv+b);
    VolT_calc = 0.5.*((pi/3).*(IDlv.^3).*L);

    % Sum of squared errors to minimized distance between measurements
    sse_val = sum((Vol_real - VolT_calc).^2);
end
