function Scatter_Vols(PIDs)
close all;
% Calculates the volumes of the LV at systole and diastole separately, with
% different parametrs for each

% Preallocate empty vectors
VolD_vec = []; % vector to store MRI max volumes
VolS_vec = []; % vector to store MRI min volumes
IDlvD_vec = []; % vector to store calculated max volumes from TTE data
IDlvS_vec = []; % vector to store calculated min volumes from TTE data
% CO_vec = []; % vector to store CO values
% HR_vec = []; % vector to store HR values

% change directory to extract patient data
cd ../Core/

for k = 1:length(PIDs)
    PID = PIDs(k); 
    data = Patient(PID); 

    % Store MRI max and min volumes to vectors
    VolD_vec(k) = data.VlvM;
    VolS_vec(k) = data.VlvmT;

    % Store TTE max and min volumes to vectors
    IDlvD_vec(k) = data.IDlvD; % inner diameter at its largest during diastole, measured in cm
    IDlvS_vec(k) = data.IDlvS; % inner diameter at its largest during systole, measured in cm

    % Store CO and HR
    % CO_vec(k) = data.CO;
    % HR_vec(k) = data.HR;
end

% back to this directory
cd ../Teichholz/

% Check that vector entries are okay
good = ~isnan(VolD_vec) & ~isnan(VolS_vec) & ~isnan(IDlvD_vec) & ~isnan(IDlvS_vec);
VolD_vec = VolD_vec(good);
VolS_vec = VolS_vec(good);
IDlvD_vec = IDlvD_vec(good);
IDlvS_vec = IDlvS_vec(good);

% Initiual guess for a and b (used for both diastolic and systolic volume)
guess = [0.0425, 0.125]; % guess and check values

% Parameter optimization for a and b 
% DO SEPARATE OPTIMIZATION FOR SYSTOLIC AND DIASTOLIC
parmsD = fminsearch(@(p) vol_errorD(p, IDlvD_vec, VolD_vec), guess);
parmsD_a = parmsD(1);
parmsD_b = parmsD(2);
disp(parmsD_a)
disp(parmsD_b)

parmsS = fminsearch(@(p) vol_errorS(p, IDlvS_vec, VolS_vec), guess);
parmsS_a = parmsS(1);
parmsS_b = parmsS(2);
disp(parmsS_a)
disp(parmsS_b)

% Calculations of volume are based on TTE data using Teichholz equation
LD = 1./(parmsD_a.*IDlvD_vec+parmsD_b); % L/D for diastole
VolTD = 0.5.*((pi/3).*(IDlvD_vec.^3).*LD); % calculated TTE volume at diastole

LS = 1./(parmsS_a.*IDlvS_vec+parmsS_b); % L/D for systole
VolTS = 0.5.*((pi/3).*(IDlvS_vec.^3).*LS); % calculated TTE volume at systole

% Linear regression
% Calculating slope and intercept for diastolic volume
coeffsD = polyfit(VolD_vec, VolTD, 1);
mD = coeffsD(1);

% Calculating slope and intercept for systolic volume
coeffsS = polyfit(VolS_vec, VolTS, 1);
mS = coeffsS(1);

% To graph y = x
xRangeD = linspace(min(VolD_vec), max(VolD_vec), 100);
xRangeS = linspace(min(VolS_vec), max(VolS_vec), 100);

%% To graph Bland-Altman Plots
avgD = (VolTD+VolD_vec)/2; 
avgS = (VolTS+VolS_vec)/2;

diffD = VolD_vec - VolTD; 
diffS = VolS_vec - VolTS; 

biasD = mean(diffD);
biasS = mean(diffS);

sdDiffD = std(diffD);
sdDiffS = std(diffS);

uLOAD = biasD + (1.96*sdDiffD);
lLOAD = biasD - (1.96*sdDiffD);
uLOAS = biasS + (1.96*sdDiffS);
lLOAS = biasS - (1.96*sdDiffS);

%% Plots
% Plotting the scatter plots
figure(1);
subplot(1,2,1); % diastolic scatter plot
    scatter(VolD_vec, VolTD, 'filled', 'MarkerFaceColor','b'); % data and calculated values
    hold on;
    y_fitD = polyval(coeffsD, VolD_vec);
    plot(VolD_vec, y_fitD, 'r', 'LineWidth', 2, 'LineStyle', '-'); % regression
    plot(xRangeD, xRangeD, 'm', 'LineWidth',2, 'LineStyle', ':'); % y = x
    hold off;
    title('MRI vs TTE Left Ventricular Volume at Diastole');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend('Data Points', 'Regression Line', 'y=x','Location', 'southeast');
    slopeText = sprintf('Slope = %.2f', mD);
    annotation('textbox', [.15 .7 .15 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
subplot(1,2,2); % systolic scatter plot
    scatter(VolS_vec, VolTS, 'filled', 'MarkerFaceColor','b');
    hold on;
    y_fitS = polyval(coeffsS, VolS_vec);
    plot(VolS_vec, y_fitS, 'r', 'LineWidth', 2, 'LineStyle', '-');
    plot(xRangeS, xRangeS, 'm', 'LineWidth',2, 'LineStyle', ':');
    hold off;
    title('MRI vs TTE Left Ventricular Volume at Systole');
    xlabel('MRI Volume (mL)');
    ylabel('TTE Volume (mL)');
    legend('Data Points', 'Regression Line', 'y=x','Location', 'southeast');
    slopeText = sprintf('Slope = %.2f', mS);
    annotation('textbox', [.60 .7 .1 .13], 'String', slopeText, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
% Bland-Altman Plots
figure(2);
subplot(1,2,1); % diastolic Bland-Altman plot
    plot(avgD, diffD, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    hold on;
    xlims = xlim;
    line(xlims, [biasD, biasD], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOAD, uLOAD], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOAD, lLOAD], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasD + (max(diffD)*0.05), sprintf('Bias: %.2f', biasD), 'Color', 'r');
    text(xlims(2)*0.75, uLOAD + (max(diffD)*0.05), sprintf('+1.96 SD: %.2f', uLOAD), 'Color', 'r');
    text(xlims(2)*0.75, lLOAD - (max(diffD)*0.05), sprintf('-1.96 SD: %.2f', lLOAD), 'Color', 'r');
    grid on;
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE (MRI - TTE)');
    title('Bland-Altman Plot for Diastolic Volume of LV');
    hold off;
subplot(1,2,2); % systolic Bland-Altman plot
    plot(avgS, diffS, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    hold on;
    xlims = xlim;
    line(xlims, [biasS, biasS], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    line(xlims, [uLOAS, uLOAS], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    line(xlims, [lLOAS, lLOAS], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    text(xlims(2)*0.75, biasS + (max(diffS)*0.05), sprintf('Bias: %.2f', biasS), 'Color', 'r');
    text(xlims(2)*0.75, uLOAS + (max(diffS)*0.05), sprintf('+1.96 SD: %.2f', uLOAS), 'Color', 'r');
    text(xlims(2)*0.75, lLOAS - (max(diffS)*0.05), sprintf('-1.96 SD: %.2f', lLOAS), 'Color', 'r');
    grid on;
    xlabel('Mean of MRI and TTE Values');
    ylabel('Difference Between MRI and TTE (MRI - TTE)');
    title('Bland-Altman Plot for Systolic Volume of LV');
    hold off;
% Box plots
figure(3);
subplot(1,2,1);
    b = boxplot(diffD);
    xlabel('Diastolic LV')
    ylabel('Difference Between MRI and TTE (MRI - TTE)')
    set(b, 'linewidth',2)
subplot(1,2,2);
    h = boxplot(diffS);
    xlabel('Systolic LV')
    set(h, 'LineWidth', 2)
end

% Nested error function for diastolic volume
function sse_D = vol_errorD(params, IDlvD, VolD_real)
    a = params(1);
    b = params(2);

    % Diastolic calculation
    LD = 1 ./ (a .* IDlvD + b);
    VolTD_calc = 0.5 .* ((pi/3) .* (IDlvD.^3) .* LD);

    % Calculate the actual regression slopes for this specific a and b
    % coeffsD = polyfit(VolD_real, VolTD_calc, 1);
    % mD = coeffsD(1); % Diastolic Slope
    
    sse_D = sum((VolD_real - VolTD_calc).^2);
end

% Nested error function for systolic volume
function sse_S = vol_errorS(params, IDlvS, VolS_real)
    a = params(1);
    b = params(2);

    % Systolic calculation
    LS = 1 ./ (a .* IDlvS + b);
    VolTS_calc = 0.5 .* ((pi/3) .* (IDlvS.^3) .* LS);

    % Calculate the actual regression slopes for this specific a and b
    % coeffsS = polyfit(VolS_real, VolTS_calc, 1);
    % mS = coeffsS(1); % Systolic Slope

    sse_S = sum((VolS_real - VolTS_calc).^2);
end