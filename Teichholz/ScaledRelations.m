function ScaledRelations(PIDs)
close all;

% Looks at the allometric relationship between height & volume and weight
% & volume. The allometric relationship can be described as
% Y = aX^b
% where Y is heart volume, X is the measurement being used (height or
% weight), a is the allometric coefficient, and b is the scaling exponent.

% Combines height and weight into Body Surface Area (BSA) to
% compare to volume. To calculate BSA, we will use Du Bois' equation:
% BSA = 0.007184*(W^0.425)*(H^0.725)

% Since we only have the inner diameter of the left ventricle, we will only
% carry out these comparisons for the left ventricle.

% Preallocate empty vectors for speed
VlvS = [];
VlvD = [];
IDlvS = [];
IDlvD = [];
hgt = [];
wgt = [];

% change directory to extract patient data
cd ../Core/

for i = 1:length(PIDs)
    PID = PIDs(i);
    data = Patient(PID);

    VlvS(i) = data.VlvmT;
    VlvD(i) = data.VlvM;
    IDlvD(i) = data.IDlvD;
    IDlvS(i) = data.IDlvS;
    hgt(i) = data.H;
    wgt(i) = data.W;
end

% change back directory
cd ../Teichholz/

% Check that vector entries are okay
good = ~isnan(VlvS) & ~isnan(VlvD) & ~isnan(IDlvD) & ~isnan(IDlvS) & ~isnan(hgt) & ~isnan(wgt);
VlvS = VlvS(good);
VlvD = VlvD(good);
IDlvD = IDlvD(good);
IDlvS = IDlvS(good);
hgt = hgt(good);
wgt = wgt(good);

%% Calculations and figures for BSA
BSA = 0.007184.*(wgt.^(0.425)).*(hgt.^(0.725));

figure(1); % BSA vs diastolic LV Volume
scatter(BSA, VlvD, 'filled', 'MarkerFaceColor', [0 0.5 0]);
hold on;
coeffs_BSAD = polyfit(BSA, VlvD, 1);
y_fit_BSA = polyval(coeffs_BSAD, BSA);
plot(BSA, y_fit_BSA, 'r', 'LineWidth',2, 'LineStyle', '-');
hold off;
title('BSA vs LV Diastolic Volume');
xlabel('BSA (m^2)');
ylabel('Volume (mL)');

figure(2); % BSA vs systolic LV volume
scatter(BSA, VlvS, 'filled', 'MarkerFaceColor', [0 0.5 0]);
hold on;
coeffs_BSAS = polyfit(BSA, VlvS, 1);
y_fit_BSAS = polyval(coeffs_BSAS, BSA);
plot(BSA, y_fit_BSAS, 'r', 'LineWidth',2, 'LineStyle','-');
hold off;
title('BSA vs LV Systolic Volume');
xlabel('BSA (m^2)');
ylabel('Volume (mL)');

%% Allometric Relationship between height and LV Volume
% Combining vectors to put diastolic and systolic volume on the same graph
height = [hgt, hgt];
Vols = [VlvS, VlvD];

log_H = log(height);
sorted_log_H  = sort(log_H);
log_V = log(Vols);

figure(3); % scatter plot of height vs LV volume
C = [ones(size(VlvS)), 2*ones(size(VlvD))];
scatter(log_H, log_V, 36, C, 'filled');
colormap([0 0 0.5; 0.5 0.5 1]);
hold on;
% Linear regression
coeffs_H = polyfit(log_H, log_V, 1);
disp('Slope and y-intercept:')
disp([coeffs_H(1), coeffs_H(2)])
y_fit_H = polyval(coeffs_H, sorted_log_H);
reg_H = plot(sorted_log_H, y_fit_H, 'r', 'LineWidth', 2, 'LineStyle', '-');
d_H = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
h_H = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
hold off;
title('Allometric Relationship Between Height and Volume');
xlabel('log(height) (cm)');
ylabel('log(volume) (mL)');
legend([d_H, h_H, reg_H], {'LV Sys', 'LV Dias', 'Linear Regression'}, 'Location', 'southeast');

% Convert:
disp('Allometric Coeff for Height and Volume:')
allo_a_H = 10^(coeffs_H(1));
disp(allo_a_H)

% Interpretations:
allo_b_H = coeffs_H(2);
if allo_b_H == 1
    disp('Isometry: Perfect scaling')
elseif allo_b_H > 1
    disp('Positive allometry')
elseif allo_b_H < 1
    disp('Negative allometry')
end

%% Allometric Relationship Between Weight and Volume
% Combining vectors to put diastolic and systolic volume on the same graph
weight = [wgt, wgt];

log_W = log(weight);
sorted_log_W  = sort(log_W);

figure(4); % scatter plot of height vs LV volume
C = [ones(size(VlvS)), 2*ones(size(VlvD))];
scatter(log_W, log_V, 36, C, 'filled');
colormap([0 0 0.5; 0.5 0.5 1]);
hold on;
% Linear regression
coeffs_W = polyfit(log_W, log_V, 1);
disp('Slope and y-intercept:')
disp([coeffs_W(1), coeffs_W(2)])
y_fit_W = polyval(coeffs_W, sorted_log_W);
reg_W = plot(sorted_log_W, y_fit_W, 'r', 'LineWidth', 2, 'LineStyle', '-');
d_W = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
h_W = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
hold off;
title('Allometric Relationship Between Weight and Volume');
xlabel('log(weight) (kg)');
ylabel('log(volume) (mL)');
legend([d_W, h_W, reg_W], {'LV Sys', 'LV Dias', 'Linear Regression'}, 'Location', 'southeast');

% Convert:
disp('Allometric Coeff for Weight and Volume:')
allo_a_W = 10^(coeffs_W(1));
disp(allo_a_W)

% Interpretations:
allo_b_W = coeffs_W(2);
if allo_b_W == 1
    disp('Isometry: Perfect scaling')
elseif allo_b_W > 1
    disp('Positive allometry')
elseif allo_b_W < 1
    disp('Negative allometry')
end


%% Allometric Relationship between BSA and Volume
% Combining vectors to put diastolic and systolic volume on the same graph
body = [BSA, BSA];

log_BSA = log(body);
sorted_log_B  = sort(log_BSA);

figure(5); % scatter plot of height vs LV volume
C = [ones(size(VlvS)), 2*ones(size(VlvD))];
scatter(log_BSA, log_V, 36, C, 'filled');
colormap([0 0 0.5; 0.5 0.5 1]);
hold on;
% Linear regression
coeffs_BSA = polyfit(log_BSA, log_V, 1);
disp('Slope and y-intercept:')
disp([coeffs_BSA(1), coeffs_BSA(2)])
y_fit_BSA = polyval(coeffs_BSA, sorted_log_B);
reg_BSA = plot(sorted_log_B, y_fit_BSA, 'r', 'LineWidth', 2, 'LineStyle', '-');
d_B = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'none');
h_B = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
hold off;
title('Allometric Relationship Between BSA and Volume');
xlabel('log(BSA) (m^2)');
ylabel('log(volume) (mL)');
legend([d_B, h_B, reg_BSA], {'LV Sys', 'LV Dias', 'Linear Regression'}, 'Location', 'southeast');

% Convert:
disp('Allometric Coeff for BSA and Volume:')
allo_a_B = 10^(coeffs_BSA(1));
disp(allo_a_B)

% Interpretations:
allo_b_B = coeffs_BSA(2);
if allo_b_B == 1
    disp('Isometry: Perfect scaling')
elseif allo_b_B > 1
    disp('Positive allometry')
elseif allo_b_B < 1
    disp('Negative allometry')
end