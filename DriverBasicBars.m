function  [patientTable, S_history] = DriverBasicBars(PID,opt,teich,k,ColorMap,FS)
%fuunctinname: DriverBasic
%Inputs patient ID (PID), flag 1 reads optimized parameters, 0 uses nominal
%parameters, teich 1 reads calculated volumes, 0 reads measured volumes
%integer k, the iteration number
%Requires: Patient.m (reads the patient files), load_global.m (generates
%nominal parameter values and initial conditions), solveModel.m (generates
%solution of the ODEs), this function requires modelBasic.m (containing the
%rhs of the ODEs).
%Function: Solves the ODEs with nominal (opt=0) or optimized (opt=1)
%parameters. Plot results in figure (1). If run for multiple iterations,
%the figure will overwrite the results

close all;

% Load patient data for patient PID
cd Core
data = Patient(PID, teich);

data.ODE_TOL  = 1e-10;              %ODE_TOL;
data.DIFF_INC = sqrt(data.ODE_TOL); %DIFF_INC;

%Load nominal parameter values (pars) & (parsG), and initial conditions for the ODEs
[pars parsG times Init] = load_global(PID, data);

data.parsG = parsG;
data.Init  = Init;
data.dt    = times.dt;
data.T     = times.T;
data.CC    = 5;

if opt == 1 % Load optimized parameters
    cd(sprintf('../MultiOptPatients/Patient_%d',PID));
    s = sprintf('Opt_Run%d.mat',k);
    load(s)
    INDMAP = data.INDMAP; 
    pars = optPars;
    cd ../../Core
else
    k = "not opt";
    ColorMap = [0 0 1];
end

assignin('base', 'pars_nominal', exp(pars));

% Solve the ODEs, return residual, least squares cost, time, MSE, States
data.computeAnalyticalSens = true;
[rout, J, tAll, MSE, States, S_history, S_time] = SolveModel(pars,data);

Ie = length(tAll);
Is = ceil(Ie - data.T*data.CC/data.dt);
tAllConv = tAll(Is:Ie);
tAllConv = tAllConv-tAllConv(1);

% Calculate cardiac output
CO = trapz(tAll,States.Qsys)/tAll(end)*60/1000;
CO = CO*ones(size(tAll));

cd ..

tend = tAllConv(end);
%plots to display P and V solutions
figure(1);hold on;
fontsize(FS+2,"points")
% Pressure vs Time
nexttile(1);
    plot(tAllConv,States.Ppa(Is:Ie),'color', ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tend];
    plot(xLimits, [data.PpaS, data.PpaS], 'k:', 'linewidth', 2); 
    % plotting range around systolic data line
        yL = data.PpaS-5; % lower line
        yU = data.PpaS+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    plot(xLimits, [data.PpaD, data.PpaD], 'k:', 'linewidth', 2); 
    % plotting range around diastolic data line
        yL = data.PpaD-5; % lower line
        yU = data.PpaD+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    ylim([0.95*min(min(States.Ppa(Is:Ie)), data.PpaD-5) 1.05*max(max(States.Ppa(Is:Ie)), data.PpaS+5)]);
    ylabel('Ppa (mmHg)');
nexttile(2)
    plot(tAllConv, States.Ppv(Is:Ie),'color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tend];
    plot(xLimits, [data.Ppv, data.Ppv], 'k:', 'linewidth', 2);
    % plotting range around data line
        yL = data.Ppv-5; % lower line
        yU = data.Ppv+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    ylim([0.95*min(min(States.Ppv(Is:Ie)), data.Ppv-5) 1.05*max(max(States.Ppv(Is:Ie)), data.Ppv+5)]);    
    ylabel('Ppv (mmHg)');
nexttile(3)
    plot(tAllConv, States.Plv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tend]; 
    ylim([0.95*min(States.Plv(Is:Ie)), 1.05*max(States.Plv(Is:Ie))]);    
    ylabel('Plv (mmHg)');
nexttile(4)
    plot(tAllConv, States.Psa(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tend];
    plot(xLimits, [data.PsaS, data.PsaS], 'k:', 'linewidth', 2);
    % plotting range around systolic data line
        yL = data.PsaS-5; % lower line
        yU = data.PsaS+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    plot(xLimits, [data.PsaD, data.PsaD], 'k:', 'linewidth', 2);
    % plotting range around diastolic data line
        yL = data.PsaD-5; % lower line
        yU = data.PsaD+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    ylim([0.95*min(min(States.Psa(Is:Ie)), data.PsaD-5) 1.05*max(max(States.Psa(Is:Ie)), data.PsaS+5)]);
    ylabel('Psa (mmHg)');
nexttile(5)
    plot(tAllConv, States.Psv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0,tend];
    plot(xLimits, [data.Psv, data.Psv], 'r:', 'linewidth', 2);
    ylim([0.95*min(min(States.Psv(Is:Ie)),data.Psv) 1.05*max(max(States.Psv(Is:Ie)),data.Psv)]);    
    ylabel('Psv (mmHg)');
nexttile(6)
    plot(tAllConv, States.Prv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tend];
    plot(xLimits, [data.PrvS, data.PrvS], 'k:', 'linewidth', 2);
    % plotting range around systolic data line
        yL = data.PrvS-5; % lower line
        yU = data.PrvS+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    plot(xLimits, [data.PrvD, data.PrvD], 'k:', 'linewidth', 2);
    % plotting range around diastolic data line
        yL = data.PrvD-5; % lower line
        yU = data.PrvD+5; % upper line
        xCord = [xLimits(1), xLimits(2), xLimits(2), xLimits(1)];
        yCord = [yL, yL, yU, yU];
        fill(xCord, yCord, [0 0.5 0], 'FaceAlpha',0.3, 'EdgeColor', 'none');
        plot(xLimits, [yL, yL], 'Color', [0 0.5 0], 'LineWidth', 2);
        plot(xLimits, [yU, yU], 'Color', [0 0.5 0], 'LineWidth', 2);
    ylim([0.95*min(min(States.Prv(Is:Ie)), data.PrvD-5) 1.05*max(max(States.Prv(Is:Ie)), data.PrvS+5)]);
    ylabel('Prv (mmHg)');
% Volume vs Time
nexttile(7)
    plot(tAllConv, States.Vpa(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    ylim([0.95*min(States.Vpa(Is:Ie)), 1.05*max(States.Vpa(Is:Ie))]);
    ylabel('Vpa (mL)');
nexttile(8)
    plot(tAllConv, States.Vpv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    ylim([0.99*min(States.Vpv(Is:Ie)), 1.01*max(States.Vpv(Is:Ie))]);
    ylabel('Vpv (mL)');
nexttile(9)
    plot(tAllConv, States.Vlv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tAllConv(end)];
    plot(xLimits, [data.VlvM, data.VlvM], 'k:', 'LineWidth', 2);
    plot(xLimits, [data.Vlvm, data.Vlvm], 'k:', 'LineWidth', 2);
    ylim([0.95*min(min(States.Vlv(Is:Ie)), data.Vlvm) 1.05*max(max(States.Vlv(Is:Ie)), data.VlvM)]);
    xlabel('time (s)');
    ylabel('Vlv (mL)');
nexttile(10)
    plot(tAllConv, States.Vsa(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    ylim([0.99*min(States.Vsa(Is:Ie)), 1.01*max(States.Vsa(Is:Ie))]);
    xlabel('time (s)');
    ylabel('Vsa (mL)');
nexttile(11)
    plot(tAllConv, States.Vsv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    ylim([0.999*min(States.Vsv(Is:Ie)), 1.001*max(States.Vsv(Is:Ie))]);
    xlabel('time (s)');
    ylabel('Vsv (mL)');
nexttile(12)
    plot(tAllConv, States.Vrv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xLimits = [0, tend];
    plot(xLimits, [data.VrvM, data.VrvM], 'k:', 'LineWidth', 2);
    plot(xLimits, [data.Vrvm, data.Vrvm], 'k:', 'LineWidth', 2);
    ylim([0.95*min(min(States.Vrv(Is:Ie)), data.Vrvm) 1.05*max(max(States.Vrv(Is:Ie)), data.VrvM)]);
    xlabel('time (s)');
    ylabel('Vrv (mL)');
% Pressure-Volume Loops
nexttile(13)
    plot(States.Vlv(Is:Ie),States.Plv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xlim([0.95*min(States.Vlv(Is:Ie)) 1.05*max(States.Vlv(Is:Ie))]);
    ylim([0.95*min(States.Plv(Is:Ie)) 1.05*max(States.Plv(Is:Ie))]);
    xlabel('Vlv (ml)');
    ylabel('Plv (mmHg)');
nexttile(14)
    plot(States.Vrv(Is:Ie),States.Prv(Is:Ie),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    xlim([0.95*min(States.Vrv(Is:Ie)) 1.05*max(States.Vrv(Is:Ie))]);
    ylim([0.95*min(States.Prv(Is:Ie)) 1.05*max(States.Prv(Is:Ie))]);
    xlabel('Vrv (ml)');
    ylabel('Prv (mmHg)');
% Cardiac Output
nexttile(15)
    nCycles = length(States.CO);
    cycIdxStart = max(1, nCycles - data.CC);
    cycRange = cycIdxStart:nCycles;
    cycTimes = (cycRange - cycIdxStart) * data.T;
    plot(cycTimes, States.CO(cycRange),'Color',ColorMap,'linewidth', 2); hold on;
    fontsize(FS,"points")
    plot(tAllConv, data.CO*ones(size(tAllConv)),'k:','linewidth',2)
    ylim([0.95*min(min(States.CO(cycRange)), data.CO), 1.05*max(max(States.CO(cycRange)), data.CO)]);
    xlabel('time (s)');
    ylabel('CO (L/min)')
    
hold(nexttile(1), 'on'); 
hModel       = plot(nan, nan, 'b-', 'LineWidth', 2);
hData        = plot(nan, nan, 'k:', 'LineWidth', 2);
hCalculated  = plot(nan, nan, 'r:', 'LineWidth', 2);
hError       = plot(nan, nan, 'Color', [0 0.5 0], 'LineWidth', 2);
lgd = legend([hModel, hData, hCalculated, hError], ...
    {'Model', 'Data', 'Calculated Value', 'Error Bars'}, ...
    'Orientation', 'horizontal','fontsize',FS+2);
lgd.Layout.Tile = 'south';

% Parameters
pars = exp(pars);

% Resistances
Rpul   = pars(1);
Rmval  = pars(2);
Raval  = pars(3);
Rsys   = pars(4);
Rtval  = pars(5);
Rpval  = pars(6);

% Compliances
Cpv = pars(7);
Csa	= pars(8);
Csv	= pars(9);
Cpa = pars(10);

% Heart parameters
ElvM = pars(11);
Elvm = pars(12);
ErvM = pars(13);
Ervm = pars(14);

TC = pars(15); 
TR = pars(16);

COModel = mean(States.CO(cycRange));

%disp('Check')
%disp([States.VlvM(end)-States.Vlvm(end),States.VrvM(end)-States.Vrvm(end),States.SV(end),data.SV,States.CO(end),data.CO]);

% Formulate table 
patientTable = table(PID, Rpul, Rmval, Raval, Rsys, Rtval, Rpval, Cpv, Csa, ... %9
    Csv, Cpa, ElvM, Elvm, ErvM, Ervm, COModel, MSE(1), MSE(2), MSE(3), MSE(4), ... %11
    MSE(5), MSE(6), MSE(7), MSE(8), MSE(9), MSE(10), MSE(11), 0,  J, ... %9 %0 instead of MSE(12) because CO is currently messed up
    'VariableNames', {'PID', 'Rpul', 'Rmval', 'Raval', 'Rsys', 'Rtval', 'Rpval', 'Cpv', 'Csa', 'Csv', 'Cpa', ... %11
                     'ElvM', 'Elvm', 'ErvM', 'Ervm', 'COModel',... %5
                     'PpaSmse', 'PpaDmse', 'PsaSmse', 'PsaDmse', 'PpvMmse', 'PrvDmse', 'PrvSmse', ... %
                     'VrvDmse', 'VrvSmse', 'VlvDmse', 'VlvSmse', 'COmse', 'J'}); % squared errors of the minimum 
end 

