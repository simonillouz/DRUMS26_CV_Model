function StrokeVol(PIDs)
% Calculating ESV based on EDV and SV 
% SV = CO/HR
% ESV = EDV - SV
close all;

% PREALLOCATED EMPTY VECTORS
VolD_vec = []; % vector to store MRI max volumes
VolS_vec = []; % vector to store MRI min volumes
IDlvD_vec = []; % vector to store max inner diameter from TTE data
CO_vec = []; % vector to store CO values
HR_vec = []; % vector to store HR values
SV_vec = []; % vector to store SV values

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

    % Store CO, HR, and SV to vectors
    CO_vec(k) = data.CO;
    HR_vec(k) = data.HR;
    SV_vec(k) = data.SV;
end

% back to this directory
cd ../Teichholz/

% CALCULATE DIASTOLIC VOLUME BASED ON TEICHHOLZ EQN, INCLUDING PARAMETER OPTIMZATION
guess = [0.0767, -0.0550];

pars = fminsearch(@(p) vol_err(p, IDlvD_vec, VolD_vec), guess);
a = pars(1);
b = pars(2);
disp('Parameters:')
disp(a)
disp(b)

L = 1./(a.*IDlvD_vec+b);
Vdias = 0.5.*((pi/3).*(IDlvD_vec.^3).*L);

% STROKE VOLUME CALCULATIONS BASED ON DIASTOLIC VOLUME
% SV = CO_vec./HR_vec;
Vsys = Vdias - SV_vec;
disp(Vsys)
pause;

% CALCULATIONS OF SLOPE FOR LINEAR REGRESSION
validD_fit = isfinite(VolD_vec) & isfinite(Vdias);
coeffsD = polyfit(VolD_vec(validD_fit), Vdias(validD_fit), 1);
mD = coeffsD(1);
mdlD = fitlm(VolD_vec(validD_fit), Vdias(validD_fit));
R2D2 = mdlD.Rsquared.Ordinary;
disp('R^2 Value for Diastolic Volumes:')
disp(R2D2)
 
validSV_fit = isfinite(VolS_vec) & isfinite(Vsys);
coeffsSV = polyfit(VolS_vec(validSV_fit), Vsys(validSV_fit), 1);
mSV = coeffsSV(1);
mdlSV = fitlm(VolS_vec(validSV_fit), Vsys(validSV_fit));
R2SV = mdlSV.Rsquared.Ordinary;
disp('R^2 Value for Systolic Volumes via SV:')
disp(R2SV)

% Y = X FOR COMPARISON
xRangeD = linspace(min(VolD_vec), max(VolD_vec), 100);
xRangeSV = linspace(min(VolS_vec), max(VolS_vec), 100);

% INFO FOR BLAND-ALTMAN PLOTS
avgD = (Vdias + VolD_vec)/2;
avgSV = (Vsys + VolS_vec)/2;
diffD = VolD_vec - Vdias;
diffSV = VolS_vec - Vsys;

biasD = mean(diffD(isfinite(diffD)));
biasSV = mean(diffSV(isfinite(diffSV)));
sdDiffD = std(diffD(isfinite(diffD)));
sdDiffSV = std(diffSV(isfinite(diffSV)));
uLOAD = biasD + (1.96*sdDiffD);
lLOAD = biasD - (1.96*sdDiffD);
uLOASV = biasSV + (1.96*sdDiffSV);
lLOASV = biasSV - (1.96*sdDiffSV);

% PLOTS
figure(1); % scatter plot
subplot(1,2,1); % diastolic volume
    scatter(VolD_vec, Vdias, 'filled', 'MarkerFaceColor','b');
    hold on;
    y_fitD = polyval(coeffsD, VolD_vec);
    plot(VolD_vec, y_fitD, 'r', 'LineWidth',2, 'LineStyle','-');
    plot(xRangeD, xRangeD, 'm', 'LineWidth', 2, 'LineStyle', ':');
    hold off;
    title('MRI vs TTE Left Ventricular Volume at Diastole');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend('Data Points', 'Regression Line', 'y=x', 'Location', 'southeast');
    slopeText = sprintf('Slope = %.2f', mD);
    annotation('textbox', [.15 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
subplot(1,2,2); % systolic volume based on SV
    validSV = isfinite(VolS_vec) & isfinite(Vsys);
    scatter(VolS_vec(validSV), Vsys(validSV), 'filled', 'MarkerFaceColor','b');
    hold on;
    y_fitSV = polyval(coeffsSV, VolS_vec(validSV));
    plot(VolS_vec(validSV), y_fitSV, 'r', 'LineWidth',2, 'LineStyle','-');
    plot(xRangeSV, xRangeSV, 'm', 'LineWidth', 2, 'LineStyle', ':');
    hold off;
    title('MRI vs SV Left Ventricular Volume at Systole');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend('Data Points', 'Regression Line', 'y=x', 'Location', 'southeast');
    slopeText = sprintf('Slope = %.2f', mSV);
    annotation('textbox', [.60 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
 
figure(2); % Bland-Altman Plot
subplot(1,2,1); % diastolic volume
    plot(avgD(isfinite(diffD)), diffD(isfinite(diffD)), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    hold on;
    xlims = xlim;
    line(xlims, [biasD, biasD], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOAD, uLOAD], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOAD, lLOAD], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasD + (max(diffD(isfinite(diffD)))*0.05), sprintf('Bias: %.2f', biasD), 'Color', 'r');
    text(xlims(2)*0.75, uLOAD + (max(diffD(isfinite(diffD)))*0.05), sprintf('+1.96 SD: %.2f', uLOAD), 'Color', 'r');
    text(xlims(2)*0.75, lLOAD - (max(diffD(isfinite(diffD)))*0.05), sprintf('-1.96 SD: %.2f', lLOAD), 'Color', 'r');
    grid on;
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE (MRI - TTE)');
    title('Diastolic Volume of LV');
    hold off;
subplot(1,2,2); % systolic volume based on SV
    validBA = isfinite(avgSV) & isfinite(diffSV);
    plot(avgSV(validBA), diffSV(validBA), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    hold on;
    xlims = xlim;
    line(xlims, [biasSV, biasSV], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOASV, uLOASV], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOASV, lLOASV], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasSV + (max(diffSV(validBA))*0.05), sprintf('Bias: %.2f', biasSV), 'Color', 'r');
    text(xlims(2)*0.75, uLOASV + (max(diffSV(validBA))*0.05), sprintf('+1.96 SD: %.2f', uLOASV), 'Color', 'r');
    text(xlims(2)*0.75, lLOASV - (max(diffSV(validBA))*0.05), sprintf('-1.96 SD: %.2f', lLOASV), 'Color', 'r');
    grid on;
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE (MRI - TTE)');
    title('Systolic Volume of LV from SV Calculation');
    hold off;
 
figure(3); % box and whisker plots
subplot(1,2,1); % diastolic volume
    bx1 = boxplot(diffD(isfinite(diffD)));
    ylabel('Difference Between MRI and TTE');
    title('LV Diastolic Volume');
    set(bx1, 'LineWidth', 2);
subplot(1,2,2); % systolic volume based on SV
    bx2 = boxplot(diffSV(isfinite(diffSV)));
    ylabel('Difference Between MRI and Volume from SV');
    title('LV Systolic Volume')
    set(bx2, 'LineWidth',2)
end

%% NESTED ERROR FUNCTION
function sse = vol_err(parms, IDlvD, VolD_real)
    best_a = parms(1);
    best_b = parms(2);

    % Diastolic volume calculation
    LD = 1 ./ (best_a .* IDlvD + best_b);
    VolTD = 0.5 .* ((pi/3) .* (IDlvD.^3) .* LD);

    % Sum of squared errors to minimized distance between measurements
    sse = sum((VolD_real - VolTD).^2);
end

