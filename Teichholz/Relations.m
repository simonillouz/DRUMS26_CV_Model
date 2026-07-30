function Relations(PIDs)
close all;

% The purpose of this code is to determine the relationship between patient
% height & ventricular volume, weight & ventricular volume, and diameter &
% ventricular volume. 

% PREALLOCATE EMPTY VECTORS
VrvD = []; % right ventricular volume at diastole (mL)
VrvS = []; % right ventricular volume at systole (mL)
VlvD = []; % left ventricular volume at diastole (mL)
VlvS = []; % left ventricular volume at systole (mL)

hgt = []; % height (cm)
wgt = []; % weight (kg)

IDlvD = []; % inner diameter of the left ventricle at diastole (cm)
IDlvS = []; % inner diamter of the left ventricle at systole (cm)

sex = []; % 1 for male, 2 for female

% Change directory to extract patient data
cd ../Core/

% LOOP TO GO THOUGH EACH OF THE PATIENTS
for i = 1:length(PIDs)
    % Exract the data from the patient file
    PID = PIDs(i);
    data = Patient(PID);

    % Store data into their corresponding vectors
    VrvD(i) = data.VrvM;
    VrvS(i) = data.VrvmT;
    VlvD(i) = data.VlvM;
    VlvS(i) = data.VlvmT;

    hgt(i) = data.H;
    wgt(i) = data.W;

    IDlvD(i) = data.IDlvD;
    IDlvS(i) = data.IDlvS;

    sex(i) = data.Sex;
end

% Change directory back
cd ../Teichholz/

% CREATE PLOTS
figure(1); % height vs volume, with indicators for sex - need to debug this plot
    hold on;
    h = gscatter(hgt, VrvD, sex, ['r','r'], ['o','x']);
    j = gscatter(hgt, VrvS, sex, ['b','b'], ['o','x']);
    l = gscatter(hgt, VlvD, sex, ['g','g'], ['o','x']);
    n = gscatter(hgt, VlvS, sex, ['m','m'], ['o','x']);
    for k = 1:length(h)
            set(h(k), 'LineWidth', 2);
            set(h(k), 'MarkerSize', 5);
            set(j(k), 'LineWidth', 2);
            set(j(k), 'MarkerSize', 5);
            set(l(k), 'LineWidth', 2);
            set(l(k), 'MarkerSize', 5);
            set(n(k), 'LineWidth', 2);
            set(n(k), 'MarkerSize', 5);
    end
    xlabel('Height (cm)');
    ylabel('Ventricular Volume (mL)');
    title('Height vs Volume by Sex');
    grid on;
    legend({'RV Diastole M','RV Diastole F','RV Systole M','RV Systole F', ...
            'LV Diastole M','LV Diastole F','LV Systole M','LV Systole F'});
    hold off;

figure(2); % weight vs volume, with indicators for sex
    hold on;
    c = gscatter(wgt, VrvD, sex, ['r','r'], ['o','x']);
    d = gscatter(wgt, VrvS, sex, ['b','b'], ['o','x']);
    e = gscatter(wgt, VlvD, sex, ['g','g'], ['o','x']);
    v = gscatter(wgt, VlvS, sex, ['m','m'], ['o','x']);
    for k = 1:length(c)
            set(c(k), 'LineWidth', 2);
            set(c(k), 'MarkerSize', 5);
            set(d(k), 'LineWidth', 2);
            set(d(k), 'MarkerSize', 5);
            set(e(k), 'LineWidth', 2);
            set(e(k), 'MarkerSize', 5);
            set(v(k), 'LineWidth', 2);
            set(v(k), 'MarkerSize', 5);
    end
    xlabel('Weight (kg)');
    ylabel('Ventricular Volume (mL)');
    title('Weight vs Volume by Sex');
    grid on;
    legend({'RV Diastole M','RV Diastole F','RV Systole M','RV Systole F', ...
            'LV Diastole M','LV Diastole F','LV Systole M','LV Systole F'});
    hold off;

[IDlvD_sort, sortIdx] = sort(IDlvD);

figure(3); % diameter vs volume, with no indicator for sex - no logarithm 
scatter(IDlvD, VlvD, 'filled', 'MarkerFaceColor','b');
coeffs_scatter = polyfit(IDlvD, VlvD, 1);
coeffs_quad = polyfit(IDlvD, VlvD, 5);
y_quad = polyval(coeffs_quad, IDlvD_sort);
y_fit_scatter = polyval(coeffs_scatter, IDlvD_sort);
hold on;
plot(IDlvD_sort, y_fit_scatter, 'r', 'LineWidth', 2, 'LineStyle', '-');
plot(IDlvD_sort, y_quad, 'm', 'LineWidth',2, 'LineStyle', '-');
hold off;
xlabel('LV Diameter (cm)');
ylabel('LV Volume (mL)');
title('Diameter vs Volume at Diastole');
legend({'Data', 'Linear Regression', 'Quadratic Regression'});

% regression
disp('Linear R^2')
mdl_lin = fitlm(IDlvD, VlvD);
R2_lin = mdl_lin.Rsquared.Ordinary;
disp(R2_lin)

disp('Quadratic R^2')
mdl_quad = fitlm(IDlvD, VlvD, 'poly8');
disp(mdl_quad.Rsquared.Ordinary)

% mdl_log = fitlm(log(IDlvD), VlvD);
% disp(mdl_log.Rsquared.Ordinary)

% mdl_loglog = fitlm(log(IDlvD), log(VlvD));
% disp(mdl_loglog.Rsquared.Ordinary)

figure(4); % diameter vs volume, with no indicator for sex - taking the log of diameters
% not as strong as a relationship 
scatter(log(IDlvD), log(VlvD), 'filled', 'MarkerFaceColor', 'b');
coeffs_log = polyfit(log(IDlvD), log(VlvD), 1);
y_fit_log = polyval(coeffs_log, log(IDlvD));
hold on;
plot(log(IDlvD), y_fit_log, 'r', 'LineWidth', 2, 'LineStyle', '-');
hold off;
xlabel('log of LV Diameter (cm)');
ylabel('log of LV Volume (mL)');
title('log of Diameter vs log of Volume at Diastole');
legend({'Data', 'Regression'});
%%
coeffs = polyfit(log(hgt), log(VlvS), 1);
y_fit = polyval(coeffs, log(hgt));

figure(5); % no clear trend between height and volume
a = gscatter(log(hgt), log(VlvS), sex, ['r','b'], ['o', 'x']);
for k = 1:length(a)
    set(a(k), 'LineWidth', 2);
    set(a(k), 'MarkerSize', 5);
end
hold on;
plot(log(hgt), y_fit, 'm', 'LineWidth', 2, 'LineStyle', '-');
xlabel('Height (cm)');
ylabel('LV Volume at Systole (mL)');
title('Height vs Volume by Sex');
legend({'LV Systole M', 'LV Systole F', 'Linear Regression'});
hold off;

figure(6); % no clear relation between weight and volume
b = gscatter(log(wgt), log(VlvS), sex, ['r', 'b'], ['o', 'x']);
for k = 1:length(b)
    set(b(k), "LineWidth", 2);
    set(b(k), 'MarkerSize', 5);
end
hold on;
coeff = polyfit(log(wgt), log(VlvS), 1);
yfit = polyval(coeff, log(wgt));
plot(log(wgt), yfit, 'm', 'LineWidth', 2, 'LineStyle', '-');
xlabel('Weight (kg)');
ylabel('LV Volume at Systole (mL)');
title("Weight vs Volume by Sex");
legend({'LV Systole M', 'LV Systole F', 'Linear Regression'});
hold off

%% Looking only at VlvS wrt weight and height for male patients 
% preallocate for speed
VlvS_M = [];
wgt_M = [];
hgt_M = [];

VlvS_W = [];
wgt_W = [];
hgt_W = [];

for k = 1:length(PIDs)
    PID = PIDs(k);
    data = Patient(PID);
    sex = data.Sex;

    if sex == 1
        VlvS_M(k) = data.Vlvm;
        wgt_M(k) = data.W;
        hgt_M(k) = data.H;
    elseif sex == 2
        VlvS_W(k) = data.Vlvm;
        wgt_W(k) = data.W;
        hgt_W(k) = data.H;
    end
end

figure(7); % height vs VlvS
subplot(1,2,1); % for men
    scatter(VlvS_M, hgt_M, 'filled', 'MarkerFaceColor', [0 0.5 0]);
    % scatter(hgt_M, VlvS_M, 'filled', 'MarkerFaceColor', [0 0.5 0]);
    tallmen = polyfit(VlvS_M, hgt_M, 1);
    fittallmen = polyval(tallmen, VlvS_M);
    hold on;
    plot(VlvS_M, fittallmen, 'm', 'LineWidth', 2, 'LineStyle', '-');
    hold off;
    xlabel('Height (cm)');
    ylabel('LV Volume at Systole (mL)');
    title('Height vs LV Volume for Male Patients');
subplot(1,2,2); % for women
    scatter(VlvS_W, hgt_W, 'filled', 'MarkerFaceColor', [0 0.5 0]);
    tallwomen = polyfit(VlvS_W, hgt_W, 1);
    fittallwomen = polyval(tallwomen, VlvS_W);
    hold on;
    plot(VlvS_W, fittallwomen, 'm', 'LineWidth', 2, 'LineStyle', '-');
    hold off;
    xlabel('Height (cm)');
    ylabel('LV Volume at Systole (mL)');
    title('Height vs LV Volume for Female Patients');

figure(8); % scatter plot of weight vs VlvS
subplot(1,2,1); % for men
    scatter(VlvS_M, wgt_M, 'filled', 'MarkerFaceColor', [0 0.5 0]);
    men_W = polyfit(VlvS_M, wgt_M, 1);
    fitmen = polyval(men_W, VlvS_M);
    hold on;
    plot(VlvS_M, fitmen, 'm', 'LineWidth', 2, 'LineStyle', '-');
    hold off;
    xlabel('Weight (kg)');
    ylabel('LV Volume at Systole (mL)');
    title('Weight vs LV Volume for Male Patients');
subplot(1,2,2); % for women
    scatter(VlvS_W, wgt_W, 'filled', 'MarkerFaceColor', [0 0.5 0]);
    women_W = polyfit(VlvS_W, wgt_W, 1);
    fitwomen = polyval(women_W, VlvS_W);
    hold on;
    plot(VlvS_W, fitwomen, 'm', 'LineWidth', 2, 'LineStyle', '-');
    hold off;
    xlabel('Weight (kg)');
    ylabel('LV Volume at Systole (mL)');
    title('Weight vs LV Volume for Female Patients');

% % Display R^2 values for these four graphs
% % Men
% disp('R^2 for Height of Men')
% % mdl_hgt = fitlm(hgt_M, VlvS_M);
% mdl_hgt = fitlm(VlvS_M, hgt_M);
% disp(mdl_hgt.Rsquared.Ordinary)
% 
% disp("R^2 for Weight of Men")
% mdl_wgt = fitlm(VlvS_M, wgt_M);
% disp(mdl_wgt.Rsquared.Ordinary)
% 
% % Women
% disp('R^2 for Height of Women')
% wdl_hgt = fitlm(VlvS_W, hgt_W);
% disp(wdl_hgt.Rsquared.Ordinary)
% 
% disp('R^2 for Weight of Women')
% wdl_wgt = fitlm(VlvS_W, wgt_W);
% disp(wdl_wgt.Rsquared.Ordinary)

