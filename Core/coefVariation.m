% Get parameters postoptimization
function coefVariation(PID,teich)
close all;

pars_names = {'Rpul', 'Rmval', 'Raval', 'Rsys', 'Rtval', 'Rpval', ... % 1-6
              'Cpv',  'Csa',   'Csv',   'Cpa', ...                    % 7-10
              'ElvM', 'Elvm',  'ErvM',  'Ervm', 'TC', 'TR'};          % 11-16


if teich == 0
    cd(sprintf('../MultiOptPatients/Patient_%d',PID));
elseif teich  == 1
    cd(sprintf('../MultiOptPatientsTeich/Patient_%d',PID));
end

load('Opt_Run0.mat');
    
NumRuns = data.numMSOpt;
tol     = 1e-3;
hi    = zeros(length(data.INDMAP),NumRuns);
low   = zeros(length(data.INDMAP),NumRuns);
parsN = zeros(length(data.INDMAP),NumRuns);
parsO = zeros(length(data.INDMAP),NumRuns);

for k = 1:NumRuns

    s = sprintf('Opt_Run%d.mat',k);
    load(s);
    
    INDMAP = data.INDMAP;
    hi(:,k)     = exp(data.hi(INDMAP));
    low(:,k)    = exp(data.low(INDMAP));

    parsN(:,k) = exp(nomPars(INDMAP));
    parsO(:,k) = exp(optPars(INDMAP));
end
cd ../../Core

fid = fopen('HiLoBound.txt', 'a');

for i = 1:length(INDMAP)
    mean_pars_opt(i) = mean(parsO(i,:));
    mean_pars_nom(i) = mean(parsN(i,:));
    sd_pars_opt(i)   = std(parsO(i,:));
    
    if abs(hi(i,1)-mean_pars_opt(i)) < tol
        fprintf(fid, 'PID:%s\n', num2str(PID));
        fprintf(fid, 'hi bound:\n');
        fprintf(fid, '%g %g %g\n', INDMAP(i), mean_pars_opt(i), hi(i,1));
    end
    
    if abs(low(i,1)-mean_pars_opt(i)) < tol
        fprintf(fid, 'PID:%s\n', num2str(PID));
        fprintf(fid, 'low bound:\n');
        fprintf(fid, '%g %g %g\n', INDMAP(i), low(i,1), mean_pars_opt(i));
    end
end

% for i = 1:length(INDMAP)
%     mean_pars_opt(i) = mean(parsO(i,:));
%     mean_pars_nom(i) = mean(parsN(i,:));
%     sd_pars_opt(i)   = std(parsO(i,:));
%     if abs(hi(i,1)-mean_pars_opt(i))<tol
%         disp(strcat('PID:',num2str(PID))
%         disp('hi bound:')
%         disp([INDMAP(i) mean_pars_opt(i) hi(i,1)])
%     end
%     if abs(low(i)-mean_pars_opt(i))<tol
%         disp(strcat('PID:',num2str(PID))
%         disp('low bound:')
%         disp([INDMAP(i) low(i,1) mean_pars_opt(i)])
%     end
% end

CoV = sd_pars_opt./mean_pars_opt;

x = [1:length(INDMAP)];
figure(1); clf;
h = plot(x, CoV, 'x', 'Markersize', 10, 'linewidth', 3);
set(gca, 'FontSize', 18);
title(sprintf('All pars for PID: %d', PID));
xticks(1:length(pars_names(INDMAP)));
set(gca, 'XTickLabel', pars_names(INDMAP), 'TickLabelInterpreter', 'latex', 'fontsize', 18);
xtickangle(45);
grid on;
ylabel('Coef Var')
xlim([1 length(INDMAP)]);
print(sprintf('../coefVarFigures/CoefVar_%d.png',PID),'-dpng');
cd ../
