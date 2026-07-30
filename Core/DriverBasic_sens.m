%--------------------------------------------------------------------------
%Computes the sensitivity matrix dy/dpars.
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set 
%pars = log(pars) in DriverBasic_sens
%Plots ranked sensitivities 
%--------------------------------------------------------------------------
function Isens = DriverBasic_sens(PID, teich)

global ODE_TOL DIFF_INC
ODE_TOL = 1e-6;
DIFF_INC = sqrt(ODE_TOL);

data = Patient(PID, teich);
[pars, parsG, times, Init]  = load_global(PID, data);
data.parsG = parsG;
data.Init  = Init;
data.dt    = times.dt;
data.T     = times.T;
data.CC    = 5;

%senseq finds the non-weighted sensitivities
sens = senseq(pars,data);

% load sens.mat
% ranked classical sensitivities
[M,N] = size(sens);
for i = 1:N
  sens_norm(i)=norm(sens(:,i),2);
end

% Divide the array by its maximum value so all points scale between 0 and 1
sens_norm_scaled = sens_norm/max(sens_norm);

% Sort the scaled 0-1 array
[Rsens,Isens] = sort(sens_norm_scaled,'descend');
% display(Isens);
% disp(Rsens) % what is this?

names = {'Rpul','Rmval','Raval','Rsys','Rtval','Rpval', ... 
         'Cpv','Csa','Csv','Cpa', ... 
         'ElvM','Elvm','ErvM','Ervm','TC','TR'};
xticks(1:16);

%Ranked sensitivities
figure('Name', 'Finite Difference Parameter Sensitivity Ranking');clf; 
% Switched from semilogy to plot, explicitly matched blue X properties
plot(1:16, Rsens, 'x', 'LineWidth', 4, 'MarkerSize', 20, 'Color', '#0072BD');
grid on;

% Lock boundaries tightly between 0 and 1
xlim([0.5, 16.5]);
xticks(1:16);
ylim([0, 1.05]);

set(gca, 'XTick', 1:16, 'XTickLabel', names(Isens), 'FontSize', 20);
xtickangle(45)

% Informative, scientifically precise titles and labels
ylabel('Relative Residual Sensitivity Norm [0-1]', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('Parameters (Ranked)', 'FontSize', 20, 'FontWeight', 'bold');
title(sprintf('Patient %d: Residual Sensitivity Rankings (Finite Difference Method)', PID), 'FontSize', 20, 'FontWeight', 'bold');

% h=semilogy(Rsens./max(Rsens),'x');
% set(h,'linewidth',4);
% set(h,'Markersize',20);
% set(gca,'Fontsize',20);
% xticks([1:16]);
% xticklabels(names(Isens));
% xtickangle(45)
% xlim([0 17])
% grid on;
% ylabel('Sensitivites');
% xlabel('Parameters');
cd ..
if ~exist('SensResultsFigs', 'dir')
    mkdir('SensResultsFigs');
end
sf20 = fullfile('SensResultsFigs', strcat('Sens', num2str(PID), '.png'));

%  sf20 = strcat('Sens',num2str(PID),'.png');
print(sf20,'-dpng');
cd Core/

% Saving sensitivity matrix to a separate folder
cd ..
folderName = 'SensResults';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = sprintf('sens%d.mat', PID);
filePath = fullfile(folderName, fileName);

save(filePath);
cd Core

% ss = strcat('sens',num2str(PID),'.mat');
% save(ss); % Saves file corresponding to dataset number.
