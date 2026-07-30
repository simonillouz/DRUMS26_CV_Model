close all
ODE_TOL  = 1e-10;
DIFF_INC = sqrt(ODE_TOL);

rng(5);
PIDs = randi([1, 343], 1, 10);

%PIDs_2 = randi([1, 343], 1, 10);

BatchTableSens = table(); % Initialize empty table for sensivity rankings spreadsheet
BatchTableCorr = table();

% 13  14  1  9  11  12  4  8  7  10  15  16  6   3  5  2
badset = [3, 111, 178, 228];
for i = 1:length(PIDs) % change to PIDs_2 if needed
    PID = PIDs(i); % change to PIDs_2 if needed
    sprintf('Patient %d', PID)
    if ~(ismember(PID, badset))
       % INDMAP = [1 5 6 8 9 10 11 12 13 14 16];  
       % INDMAP = [1:12 13 14 15]; 
       Isens = DriverBasic_sens(PID, 0) % needs sens matrix before it can output 
       
       % Convert to a one-row table so it can vertcat with BatchTableSens
       IsensTable = array2table(Isens(:)', 'VariableNames', strcat('Param', string(1:numel(Isens))));
       PIDTable = table(PID, 'VariableNames', {'PID'});
       rowTable = [PIDTable, IsensTable];
       BatchTableSens = [BatchTableSens; rowTable];
       
       % Both of the lists below work, but appears necessary to remove 7 parameters
       %INDMAP = [1 4 8:14]; 
       % INDMAP = [4 6 8:14];
       % INDMAP = [1 4 8 10:14]; % - does not work, still have a list of correlated parameters
        INDMAP = [1 3:14];

       folderName = 'SensResults'; % Must match the folder name in DriverBasic_opt
       cd ..
       fileName = sprintf('sens%d.mat', PID);
       filePath = fullfile(folderName, fileName);
        
       % Check if the file exists before loading to prevent crashes
       if exist(filePath, 'file')
           load(filePath);
       else
           error('Sensivity matrix not found for Patient %d at: %s', PID, filePath);
       end
       
       cd Core/

       sens = sens(:,INDMAP)
       pause
        [m,n] = size(sens);
        A  = sens'*sens;
        Ai = inv(A);
        disp('condition number of A = transpose(S)S and S');
        disp([ cond(A) cond(sens)] );
        [a,b] = size(Ai);
        r = zeros(a,b); % reset every iteration
        for ii = 1:a
            for jj = 1:b
                r(ii,jj) = Ai(ii,jj)/sqrt(Ai(ii,ii)*Ai(jj,jj)); % covariance matrix
            end
        end
        rn = triu(r,1); % extract upper triangular part of the matrix
        [rowIdx, colIdx] = find(abs(rn) > 0.95);
        disp('correlated parameters');
        for k = 1:length(rowIdx)
            p1 = INDMAP(rowIdx(k));
            p2 = INDMAP(colIdx(k));
            corrVal = rn(rowIdx(k), colIdx(k));
            disp([p1, p2, corrVal]);
            corrRow = table(PID, p1, p2, corrVal, ...
                'VariableNames', {'PID', 'Param1', 'Param2', 'Correlation'});
            BatchTableCorr = [BatchTableCorr; corrRow];
        end
       % rn = triu(r,1); % extract upper triangular part of the matrix
       % [i,j] = find(abs(rn) > 0.95);
       % 
       % disp('correlated parameters');
       % for k = 1:length(i)
       %      disp([INDMAP(i(k)),INDMAP(j(k)),rn(i(k),j(k))]);
       % end
    end
       
        %     %% SVD
        %     %check eigen values for singularity
        %     [U,Si,V] = svd(sens, 0);  %SVD decomposition
        %     svals = diag(Si);      %vector with singular values
        %     svals./svals(1);       %normalize vector with singular values and print to screen
        %     numopt = max(find(svals./svals(1) > 10*DIFF_INC)); %Why normalize? in pope,  dont sclae
        % 
        %     % number of parameters that can be estimated, note given that
        %     % DIFF_INC = sqrt(ODE_TOL), then sensitivities are accuratly computed
        %     % (using finite differencees) to order DIFF_INC, in fact they are
        %     % (DIFF_INC [big O], i.e. sensitivites accurate to some constant
        %     % multiplied by DIFF_INC, we set the bound at 10*DIFF_INC
        % 
        %     [q,r,imap] = qr(V(:,1:numopt)',0); % Perform the QR factorization to get the permutation matrix
        %     disp('subset:');
        %     sort(INDMAP(imap(1:numopt)))       % The subset is found by taking out imap parameters out of INDMAP_SVD
        %     length(INDMAP(imap(1:numopt)))
        % end 

        %check eigen values for singularity
       [U,Si,V] = svd(sens, 0);  %SVD decomposition
       svals = diag(Si);      %vector with singular values
       svals./svals(1);       %normalize vector with singular values and print to screen
       numopt = max(find(svals./svals(1) > 10*DIFF_INC)); %Why normalize? in pope,  dont sclae
    
       [q,r,imap] = qr(V(:,1:numopt)',0); % Perform the QR factorization to get the permutation matrix
       disp('subset:'); % what is a subset of
       sort(INDMAP(imap(1:numopt)))       % The subset is found by taking out imap parameters out of INDMAP_SVD
       length(INDMAP(imap(1:numopt)))

       % pwd
       % cd Core/
end

%writetable(BatchTableSens, 'Patient Sensitivities.xlsx');
%writetable(BatchTableCorr, 'Patient Correlations_1.xlsx');
