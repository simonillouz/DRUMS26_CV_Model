function MinVol(PIDs)
close all;

% Purpose of this function is to compare methods of calculating SV
% Function takes in the patient IDs and calcuates SV in the two following ways:
% (1) SV = CO/HR
% (2) SV = max - min
% The function then compares these two methods to check for inconsistencies

% PREALLOCATE EMPTY VECTORS 
VrvD_vec = []; % vector to store MRI max volumes - right heart
VrvS_vec = []; % vector to store MRI min volumes - right heart
VlvD_vec = []; % vector to store MRI max volumes - left heart
VlvS_vec = []; % vector to store MRI min volumes - left heart

IDlvD_vec = []; % vector to store max diameter from TTE data - left heart only
IDlvS_vec = []; % vector to store min diameter from TTE data - left heart only

CO_vec = []; % vector to store CO values
HR_vec = []; % vector to store HR values

% Change directory to extract patient data
cd ../Core/

for k = 1:length(PIDs)
    PID = PIDs(k);
    data = Patient(PID);

    % Store data into their corresponding vectors
    VrvD_vec(k) = data.VrvM;
    VrvS_vec(k) = data.VrvmT;
    VlvD_vec(k) = data.VlvM;
    VlvS_vec(k) = data.VlvmT;

    IDlvD_vec(k) = data.IDlvD;
    IDlvS_vec(k) = data.IDlvS;

    CO_vec(k) = data.CO;
    HR_vec(k) = data.HR;
end

% Change directory back
cd ../Teichholz/

% CALCULATIONS FOR STROKE VOLUME 
SV_Method1 = CO_vec./HR_vec;

%% COMPARING MEASURED AND CALCULATED VOLUME FROM STROKE VOLUME
% CALCULATIONS FOR LV AND RV VOLUME - USE SV_METHOD1
Vrvm = VrvD_vec - SV_Method1;
Vlvm = VlvD_vec - SV_Method1;

Vrvm_sort = sort(Vrvm);
Vlvm_sort = sort(Vlvm);

% INFORMATION FOR PLOTS
validRV_fit = isfinite(Vrvm) & isfinite(VrvS_vec);
coeffVrvm_lin = polyfit(Vrvm(validRV_fit), VrvS_vec(validRV_fit), 1);
coeffVrvm_quad = polyfit(Vrvm(validRV_fit), VrvS_vec(validRV_fit), 2);
coeffVrvm_log = polyfit(log(Vrvm(validRV_fit)), VrvS_vec(validRV_fit), 1);

validLV_fit = isfinite(Vlvm) & isfinite(VlvS_vec);
coeffVlvm_lin = polyfit(Vlvm(validLV_fit), VlvS_vec(validLV_fit), 1);
coeffVlvm_quad = polyfit(Vlvm(validLV_fit), VlvS_vec(validLV_fit), 2);
coeffVlvm_log = polyfit(log(Vlvm(validLV_fit)), VlvS_vec(validLV_fit), 1);

% Information for Bland-Altman Plots
% Left Ventricle
avgLV = (Vlvm + VlvS_vec)/2;
diffLV = Vlvm - VlvS_vec;
biasLV = mean(diffLV(isfinite(diffLV)));
sdDiffLV = std(diffLV(isfinite(diffLV)));
uLOALV = biasLV + (1.96*sdDiffLV);
lLOALV = biasLV - (1.96*sdDiffLV);

% Right Ventricle
avgRV = (Vrvm + VrvS_vec)/2;
diffRV = Vrvm - VrvS_vec;
biasRV = mean(diffRV(isfinite(diffRV)));
sdDiffRV = std(diffRV(isfinite(diffRV)));
uLOARV = biasRV + (1.96*sdDiffRV);
lLOARV = biasRV - (1.96*sdDiffRV);

% Calculate Coefficient of Determination
disp('Coefficient of Determinations for Minimum RV Volume')
mdlVrvm_lin = fitlm(Vrvm, VrvS_vec);
R2_Vrvm_lin = mdlVrvm_lin.Rsquared.Ordinary;
disp('Coefficient of Determination for Linear Fit:')
disp(R2_Vrvm_lin)

mdlVrvm_quad = fitlm(Vrvm, VrvS_vec, 'purequadratic');
R2_Vrvm_quad = mdlVrvm_quad.Rsquared.Ordinary;
disp('Coefficient of Determination for Quadratic Fit:')
disp(R2_Vrvm_quad)

mdlVrvm_log = fitlm(log(Vrvm), VrvS_vec);
R2_Vrvm_log = mdlVrvm_log.Rsquared.Ordinary;
disp('Coefficient of Determination for Logarithmic Fit:')
disp(R2_Vrvm_log)

if R2_Vrvm_lin > R2_Vrvm_quad && R2_Vrvm_lin > R2_Vrvm_log
    disp('Linear fit has the greatest R^2 value')
elseif R2_Vrvm_quad > R2_Vrvm_log && R2_Vrvm_quad > R2_Vrvm_lin
    disp('Quadratic fit has the greastest R^2 value')
else
    disp('Logarithmic fit has the greatest R^2 value')
end

disp('Coefficient of Determinations for Minimum LV Volume')
mdlVlvm_lin = fitlm(Vlvm, VlvS_vec);
R2_Vlvm_lin = mdlVlvm_lin.Rsquared.Ordinary;
disp('Coefficient of Determination for Linear Fit:')
disp(R2_Vlvm_lin)

mdlVlvm_quad = fitlm(Vlvm, VlvS_vec, 'purequadratic');
R2_Vlvm_quad = mdlVlvm_quad.Rsquared.Ordinary;
disp('Coefficient of Determination for Quadratic Fit:')
disp(R2_Vlvm_quad)

mdlVlvm_log = fitlm(log(Vlvm), VlvS_vec);
R2_Vlvm_log = mdlVlvm_log.Rsquared.Ordinary;
disp('Coefficient of Determination for Logarithmic Fit:')
disp(R2_Vlvm_log)

if R2_Vlvm_lin > R2_Vlvm_quad && R2_Vlvm_lin > R2_Vlvm_log
    disp('Linear fit has the greatest R^2 value')
elseif R2_Vlvm_quad > R2_Vlvm_log && R2_Vlvm_quad > R2_Vlvm_lin
    disp('Quadratic fit has the greastest R^2 value')
else
    disp('Logarithmic fit has the greatest R^2 value')
end

% PLOTS
figure(1); % scatter plot for LV volume
    scatter(Vlvm, VlvS_vec, 'filled', 'MarkerFaceColor', [0 0 0.5]);
    hold on;
    y_fitLV_lin = polyval(coeffVlvm_lin, Vlvm_sort);
    y_fitLV_quad = polyval(coeffVlvm_quad, Vlvm_sort);
    y_fitLV_log = polyval(coeffVlvm_log, log(Vlvm_sort));
    regLV_lin = plot(Vlvm_sort, y_fitLV_lin, 'r', 'LineWidth',2, 'LineStyle', '-');
    regLV_quad = plot(Vlvm_sort, y_fitLV_quad, 'b', 'LineWidth', 2, 'LineStyle', '-');
    regLV_log = plot(Vlvm_sort, y_fitLV_log, 'Color', [0 0.5 0], 'LineWidth', 2, 'LineStyle', '-');
    hold off;
    title('LV Volume From SV');
    xlabel('Calculated LV Volume (mL)');
    ylabel('Measured LV Volume (mL)');
    legend([regLV_lin, regLV_quad, regLV_log], ...
        {'Linear Regression', 'Quadratic Regression', 'Logarithmic Regression'}, ...
        'Location', 'southeast');
figure(2); % scatter plot for RV volume
    scatter(Vrvm, VrvS_vec, 'filled', 'MarkerFaceColor', [0 0 0.5]);
    hold on;
    y_fitRV_lin = polyval(coeffVrvm_lin, Vrvm_sort);
    y_fitRV_quad = polyval(coeffVrvm_quad, Vrvm_sort);
    y_fitRV_log = polyval(coeffVrvm_log, log(Vrvm_sort));
    regRV_lin = plot(Vrvm_sort, y_fitRV_lin, 'r', 'LineWidth',2, 'LineStyle', '-');
    regRV_quad = plot(Vrvm_sort, y_fitRV_quad, 'b', 'LineWidth', 2, 'LineStyle', '-');
    regRV_log = plot(Vrvm_sort, y_fitRV_log, 'Color', [0 0.5 0], 'LineWidth', 2, 'LineStyle', '-');
    hold off;
    title('RV Volume From SV');
    xlabel('Calculated RV Volume (mL)');
    ylabel('Measured RV Volume (mL)');
    legend([regRV_lin, regRV_quad, regRV_log], ...
        {'Linear Regression', 'Quadratic Regression', 'Logarithmic Regression'}, ...
        'Location', 'southeast');
figure(3); % LV volume
subplot(1,2,1); % bland-altman plot
    scatter(avgLV(isfinite(diffLV)), diffLV(isfinite(diffLV)), 64, 'k', 'o', 'MarkerFaceColor', 'b');
    hold on;
    xlims = xlim;
    line(xlims, [biasLV, biasLV], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOALV, uLOALV], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOALV, lLOALV], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasLV + (max(diffLV(isfinite(diffLV)))*0.05), sprintf('Bias: %.2f', biasLV), 'Color', 'r');
    text(xlims(2)*0.75, uLOALV + (max(diffLV(isfinite(diffLV)))*0.05), sprintf('+1.96 SD: %.2f', uLOALV), 'Color', 'r');
    text(xlims(2)*0.75, lLOALV - (max(diffLV(isfinite(diffLV)))*0.05), sprintf('-1.96 SD: %.2f', lLOALV), 'Color', 'r');
    grid on;
    xlabel('Mean of Calculated and Measured Volume');
    ylabel('Difference Between Calculated and Measured Volume');
    title('Left Ventricular Minimum Volume');
    hold off;
subplot(1,2,2); % box-and-whiskers plot
    bLV = boxplot(diffLV(isfinite(diffRV)));
    title('Difference Between Calculated and Measured LV Volume');
    set(bLV, 'LineWidth', 2);
figure(4); % RV Volume
subplot(1,2,1); % bland-altman plot
    scatter(avgRV(isfinite(diffRV)), diffRV(isfinite(diffRV)), 64, 'k', 'o', 'MarkerFaceColor', 'b');
    hold on;
    xlims = xlim;
    line(xlims, [biasRV, biasRV], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOARV, uLOARV], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOARV, lLOARV], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasRV + (max(diffRV(isfinite(diffRV)))*0.05), sprintf('Bias: %.2f', biasRV), 'Color', 'r');
    text(xlims(2)*0.75, uLOARV + (max(diffRV(isfinite(diffRV)))*0.05), sprintf('+1.96 SD: %.2f', uLOARV), 'Color', 'r');
    text(xlims(2)*0.75, lLOARV - (max(diffRV(isfinite(diffRV)))*0.05), sprintf('-1.96 SD: %.2f', lLOARV), 'Color', 'r');
    grid on;
    xlabel('Mean of Calculated and Measured Volume');
    ylabel('Difference Between Calculated and Measured Volume');
    title('Right Ventricular Minimum Volume');
    hold off;
subplot(1,2,2); % box-and-whiskers plot
    bRV = boxplot(diffRV(isfinite(diffRV)));
    title('Difference Between Calculated and Measured RV Volume');
    set(bRV, 'LineWidth', 2);

