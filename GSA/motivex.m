%% Nice Test Example
clear; clc; close all;

% --- Parameter Setup ---
N_inputs = 10;  
N_samples = 5000;  

% --- Sobol Computation ---
a = zeros(N_inputs,1);
ST_uniform = zeros(N_inputs,1);
for i = 1:N_inputs
    a(i) = 1/i; % a_i coefficients
    ST_uniform(i) = 1 - (pi^2/6 - a(i))./(pi^2/6); % total sobol' indices
end

% --- Function Setups ---
generate_X = @(N, p) randn(N, p);
evaluate_Y = @(X) motivatingexample(X);

%% --- RUN TOTAL HSIC COMPUTATIONS ---
tic;
T_hsic_tuned = get_total_HSIC_indices(generate_X, evaluate_Y, N_samples, N_inputs, 1);
exec_time_tuned = toc;

tic;
T_hsic_med = get_total_HSIC_indices(generate_X, evaluate_Y, N_samples, N_inputs, 0);
exec_time_med = toc;


%% --- PLOT AND SAVE RESULTS ---

figure('Name', 'Total Sobol Indices');
bar(1:N_inputs, ST_uniform, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'TickLabelInterpreter', 'latex');
xlim([0.5 (N_inputs +0.5)]);
xticks(1:1:N_inputs);
set(gca, 'fontsize', 18);
xlabel('Input $X_d$', 'Interpreter', 'latex');
ylabel('Total Sobol Index', 'Interpreter', 'latex');
%title('Total Sobol Indices')
ylim([0, max(ST_uniform)*1.1]);
exportgraphics(gcf,'motivex10sobol.eps','ContentType','vector')
%set(gca, 'YScale', 'log');

figure('Name', 'Total HSIC Indices Tuned Bandwidth');
bar(1:N_inputs, T_hsic_tuned, 'FaceColor', [0.2 0.6 0.8]);
xlim([0.5 (N_inputs +0.5)]);
xticks(1:1:N_inputs);
set(gca, 'fontsize', 18);
set(gca, 'TickLabelInterpreter', 'latex');
ylabel('Total HSIC Sensitivity Index', 'Interpreter', 'latex');
xlabel('Input $X_d$', 'Interpreter', 'latex');
%title('Total HSIC Sensitivity Indices')
%subtitle('Tuned Bandwidth Parameters')
exportgraphics(gcf,'motivex10tuned.eps','ContentType','vector')
%set(gca, 'YScale', 'log')
fname = ['motivex_tuned' num2str(N_samples) '.mat'];
save(fname, 'N_samples', 'T_hsic_tuned', 'N_inputs', 'exec_time_tuned');
fprintf('Analysis complete in %.2f seconds. Data saved to %s\n', exec_time_tuned, fname);

figure('Name', 'Total HSIC Indices Median Bandwidth');
bar(1:N_inputs, T_hsic_med, 'FaceColor', [0.2 0.6 0.8]);
xlim([0.5 (N_inputs +0.5)]);
xticks(1:1:N_inputs);
set(gca, 'fontsize', 18);
set(gca, 'TickLabelInterpreter', 'latex');
ylabel('Total HSIC Sensitivity Index', 'Interpreter', 'latex');
xlabel('Input $X_d$', 'Interpreter', 'latex');
%title('Total HSIC Sensitivity Indices')
%subtitle('Median Heuristic Bandwidth Parameters')
exportgraphics(gcf,'motivex10med.eps','ContentType','vector')
%set(gca, 'YScale', 'log')
fname = ['motivex_med' num2str(N_samples) '.mat'];
save(fname, 'N_samples', 'T_hsic_med', 'N_inputs', 'exec_time_med');
fprintf('Analysis complete in %.2f seconds. Data saved to %s\n', exec_time_med, fname);

function Y = motivatingexample(X)
[n,p] = size(X);

a = zeros(p:1);
for i =1:p
    a(i) = 1/i;
end

Y = zeros(n,1);
for i = 1:n
    Xi = X(i,:)';
    Y(i) = a*Xi ;
end
end
