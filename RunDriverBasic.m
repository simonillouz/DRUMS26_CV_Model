function RunDriverBasic(PID,teich)
%Name: RunDriverBasic
%Input patient ID (PID);
%Teich = 1 to use teichholz volume, 0 for MRI volume
%Requires: Estimated parameters placed in the folder MultiOptPatients/
%MultiOptPatientsTeich
%Function: Reads the file with estimated parameters for each iteration

close all;
if teich == 0
    cd(sprintf('MultiOptPatients/Patient_%d',PID));
elseif teich ==1
    cd(sprintf('MultiOptPatientsTeich/Patient_%d',PID));
end
NumRuns = numel(dir('*.mat')); %1
cd ../..
ramp  = linspace(0, 0.75, NumRuns+1)';
ramp0 = linspace(0, 0, NumRuns+1)';
grayMap = [ramp0, ramp0, ramp];
figure(1);hold on;
t=tiledlayout(4,4);
FS = 12;
disp(NumRuns)
pause
for i = 1:NumRuns
    DriverBasic(PID,1,0,i-1,grayMap(i,:),FS);
    %DriverBasic(PID,1,1,i-1,grayMap(i,:),FS); % for teichholz
end
figure(1);
for i = 1:12
    ax = nexttile(i);
    xlim(ax,[0 3.5]);
end
i = 15;
ax = nexttile(i);
xlim(ax,[0 3.5]);

s = sprintf('Pressure, Volume, and PV Loops for Patient %d', PID);
sgtitle(s,'fontsize',FS+2,'fontweight','bold');

if ~exist('PatientFiguresOpt', 'dir')
    mkdir('PatientFiguresOpt');
end

print(figure(1), sprintf('PatientFiguresOpt/Patient_%d.png', PID), '-dpng', '-r300');