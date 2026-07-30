function Cluster(file_name, sheet_name, Norm_Flag)

% The Cluster function clusters patients based on their parameters via 
% hierarchical and k means clustering methods

% Input the name of your file and specify the sheet you want data read from
% Input as a string ie, 'Patient Optimized Parameters','Sheet1'

% Norm_Flag determines what type of data normalization is conducted 
% 0 for no normalization, 1 for normalization by the mean, and 2 for normalization by the standard deviation

% Goes back into directory to find spreadsheet
cd ..

% This function applies both hierarchical and k means clustering to your
% input data set

% Input your data set as a string, EX: Cluster('Optimized parameters')

data_table = readtable(file_name,'Sheet',sheet_name,'VariableNamingRule', 'preserve', 'NumHeaderLines', 1, 'ReadVariableNames', true);      

% Extracts list of patient IDs from 1st row excel header 

PatIDs = fillmissing(string(readcell(file_name, 'Range', '1:1')), 'previous');  % reads first excel row into strings filling empty cells with previous value
patient_strings = unique(PatIDs(contains(PatIDs, 'Patient')), 'stable'); % Filters to extract PIDs and preservers order
PIDs = str2double(regexp(patient_strings, '\d+', 'match', 'once')); % finds the digits and turns the string into a PID vector
PIDs = PIDs';
 
% Getting corresponding PID Ejection Fractions (EF) from patient data 
cd ../'RunModel 6'/Core/
patientdata_table = readtable('AllPatsMRI_2_REU.xlsx','Sheet','Sheet1','VariableNamingRule','preserve');
[~, PID_indices] = ismember(PIDs, patientdata_table.Pat_Num);
EF = patientdata_table.EF(PID_indices);
cd ..
cd ../'RunModel 6'/Clustering/

% Preparing data for normalization

mean_column_flags = contains(data_table.Properties.VariableNames, 'Mean');
mean_matrix_raw = table2array(data_table(:,mean_column_flags));
A_Optp = mean_matrix_raw';
[NumPats_Optp, Num_Optp] = size(A_Optp);

% Taking the mean of each optimized parameter and centering the matrix
AMean_Optp = mean(A_Optp,1);                                % Measure mean across all patients
ACent_Optp = zeros(NumPats_Optp, Num_Optp);                 % Preallocate
    for i = 1:NumPats_Optp
        for j = 1:Num_Optp
            ACent_Optp(i,j) = A_Optp(i,j) - AMean_Optp(j);      % Centering the data matrix    
        end
    end

% Now normalizing (or not) depending on normalization flag
% Norm_Flag = 0 - no normalization 
%           = 1 - normalize by mean
%           = 2 - normalize by standard deviation
AStd_Optp = std(A_Optp);                                    % Measure std dev across all pats
ANorm_Optp = zeros(NumPats_Optp, Num_Optp);                 % Preallocate
if (Norm_Flag == 0)                                         % No normailzation
    ANorm_Optp = ACent_Optp;
elseif (Norm_Flag == 1)                                     % Normalizing by the mean
    for i = 1:NumPats_Optp
        for j = 1:Num_Optp
            ANorm_Optp(i,j) = ...
                ACent_Optp(i,j) / AMean_Optp(j);
        end
    end
elseif (Norm_Flag == 2)                                                        % Normalizing by standard deviation
    for i = 1:NumPats_Optp
        for j = 1:Num_Optp
            ANorm_Optp(i,j) = ...
                ACent_Optp(i,j) / AStd_Optp(j);
        end
    end
end

%Performing PCA
%Note that columns of AV contains the principal component values for each patient (sometimes called the score)

% Running the SVD with rows as each patient and columns as each optim param
[U_Optp,S_Optp,V_Optp] = svd(ANorm_Optp);                   % Singular value decomposition
% In this case since A = U*S*V', U will represent the rotation in the patient 
%  dimension and V' will represent the rotation in the optim param dimension.
%  Therefore our PCA score will be U*S and our PCA loadings are columns of V
AV_Optp = U_Optp * S_Optp;
% Now check to see the total variance is the same before and after the PCA
sigma_Optp = diag(S_Optp);                                  % Singular values
TotVar_Optp = norm(ANorm_Optp,'fro')^2;                     % Total variance of optim params
rho_Optp = norm(sigma_Optp)^2;  

%Perform hierarchical clustering and make a dendrogram

% Creating hierarchical tree where we use the Ward metric when joining clusters which 
%  considers the increase in the within cluster sum of squares. This is just the squared 
%  distance between each element in the cluster and it's centroid
NumHCClust_Optp = 2;
HCLink_Optp = linkage(ANorm_Optp,'ward');                   % Making the tree
HCClust_Optp = cluster(HCLink_Optp, ...                     % Cluster number to display
'Maxclust',NumHCClust_Optp); 
HCCutOff_Optp = median([HCLink_Optp(end-1,3) HCLink_Optp(end,3)]);

%Performing kmeans
% Running KMeans clustering on normalized optimized parameters

KMnsClust_Num = 2;

[KMnsClust_Optp,KMnsCntrd_Optp, ...                 % KMeans clust numbers
KMnsSumDst_Optp] = ...                          %  centroids, sum of dist
kmeans(ANorm_Optp,KMnsClust_Num, ...            %  Norm data, # of clusts 
'Distance','cityblock', ...                     %  Use L1 distance
'Display','final', ...                          %  Display results
'Replicates',20);                               %  Run 20 times  

% Cluster based on Ejection Fraction (EF)

EFclass = zeros(size(EF));
EFclass(EF < 0.5) = 1;
EFclass(EF >= 0.5) = 2;

% Getting colors 
cmap = parula(5);           % sample 5 points across parula for good spread
c1 = cmap(1,:);         % dark blue
c2 = cmap(2,:);         % blue / green
c3 = cmap(3,:);         % teal / green
c4 = cmap(4,:);         % yellow / green
c5 = cmap(5,:);         % bright yellow

% Set up single tabed figure window

MainFig = figure('Name','Clustering Results','Position',[100 100 1100 750]);
tg = uitabgroup(MainFig, 'Position',[0 0 1 1]);

% Tab 1: Dendrogram
tab1 = uitab(tg, 'Title', 'Dendrogram');
ax1 = axes('Parent', tab1);

axes(ax1)   % make ax1 the "current axes" so dendrogram draws into it
HCOptp_Dend = dendrogram(HCLink_Optp,0,'ColorThreshold', ...
    HCCutOff_Optp,'Labels',string(1:NumPats_Optp));
title(ax1, 'Hierarchical Dendrogram')

% Tab 2: HC PCA scatter
tab2 = uitab(tg, 'Title', 'PCA - Hierarchical');
ax2 = axes('Parent', tab2);
gscatter(ax2, AV_Optp(:,1), AV_Optp(:,2), HCClust_Optp, [c1; c3], '.', 35)
set(ax2, 'XTick', [], 'YTick', [], 'FontSize', 16)
xlabel(ax2, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp), 'FontSize', 18)
ylabel(ax2, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp), 'FontSize', 18)
title(ax2, 'PCA of Optimized Parameters, Clustered via Hierarchical', 'FontSize', 18)
legend(ax2, 'Cluster 1', 'Cluster 2', 'FontSize', 14)

%Tab 3: KMeans PCA scatter
tab3 = uitab(tg, 'Title', 'PCA - KMeans');
ax3 = axes('Parent', tab3);
gscatter(ax3, AV_Optp(:,1),AV_Optp(:,2),KMnsClust_Optp,[c1; c3],'.',35)
set(ax3, 'XTick', [], 'YTick', [],'FontSize',16)
xlabel(ax3, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp),'FontSize', 18)
ylabel(ax3, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp),'FontSize', 18)
title(ax3, 'PCA of Optimized Parameters, Clustered via k Means', 'FontSize', 18)
legend(ax3, 'Cluster 1', 'Cluster 2','FontSize',14)

%Tab 4: Silhouette
tab4 = uitab(tg, 'Title', 'Silhouette');
ax4 = axes('Parent', tab4);
[SilhKMns_Optp,SilhFigHndl_Optp] = silhouette(ANorm_Optp,KMnsClust_Optp,'cityblock');         
xlabel(ax4, 'Silhouette Value')
ylabel(ax4, 'Cluster')

%Tab 5: EF class scatter
tab5 = uitab(tg, 'Title', 'PCA - EF Class');
ax5 = axes('Parent', tab5);
gscatter(ax5, AV_Optp(:,1), AV_Optp(:,2), EFclass, [c1; c3], ".", 35);
set(ax5, 'XTick', [], 'YTick', [],'FontSize',16)
xlabel(ax5, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp),'FontSize',18)
ylabel(ax5, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp),'FontSize',18)
title(ax5, 'PCA of Optimized Parameters, Colored by Clinical EF Class','FontSize',18)
legend(ax5, 'HFrEF (EF<50%)','HFpEF (EF\geq50%)','FontSize',14)

%Tab 7: Agreement scatter
tab7 = uitab(tg, 'Title', 'Agreement');
ax7 = axes('Parent', tab7);

ThreeGroup = zeros(size(HCClust_Optp));
ThreeGroup(HCClust_Optp == 1 & KMnsClust_Optp == 1) = 1;
ThreeGroup(HCClust_Optp == 2 & KMnsClust_Optp == 2) = 2;
ThreeGroup(HCClust_Optp ~= KMnsClust_Optp) = 3;

gscatter(ax7, AV_Optp(:,1), AV_Optp(:,2), ThreeGroup, [c1; c3; c4], '.', 35)
set(ax7, 'XTick', [], 'YTick', [],'FontSize',16)
xlabel(ax7, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp),'FontSize', 18)
ylabel(ax7, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp), 'FontSize', 18)
title(ax7, 'Agreement Between Hierarchical and K Means Clustering','FontSize', 18)
legend(ax7, 'Agreed Cluster 1','Agreed Cluster 2','Disagree','FontSize', 14)

% Clustering Summary Table

% Getting PIDs of each cluster
HCClust1_PIDs = PIDs(HCClust_Optp == 1);
HCClust2_PIDs = PIDs(HCClust_Optp == 2);
KMClust1_PIDs = PIDs(KMnsClust_Optp == 1);
KMClust2_PIDs = PIDs(KMnsClust_Optp == 2);
EFHCClust_1 = EF(HCClust_Optp == 1);
EFHCClust_2 = EF(HCClust_Optp == 2);
EFKMnsClust_1 = EF(KMnsClust_Optp == 1);
EFKMnsClust_2 = EF(KMnsClust_Optp == 2);
MeanEF_HC1 = mean(EFHCClust_1);
MeanEF_HC2 = mean(EFHCClust_2);
MeanEF_KM1 = mean(EFKMnsClust_1);
MeanEF_KM2 = mean(EFKMnsClust_2);

% Creating Table of PIDs, Cluster, and EF - Claude

%Hierarchical: Cluster 1 and Cluster 2 side by side

HC1_PID = HCClust1_PIDs;
HC1_EF  = EFHCClust_1;
HC2_PID = HCClust2_PIDs;
HC2_EF  = EFHCClust_2;

nHC = max(length(HC1_PID), length(HC2_PID));

HC1_PID_pad = [HC1_PID; nan(nHC-length(HC1_PID),1)];
HC1_EF_pad  = [HC1_EF;  nan(nHC-length(HC1_EF),1)];
HC2_PID_pad = [HC2_PID; nan(nHC-length(HC2_PID),1)];
HC2_EF_pad  = [HC2_EF;  nan(nHC-length(HC2_EF),1)];

HC_Table = table(HC1_PID_pad, HC1_EF_pad, HC2_PID_pad, HC2_EF_pad, ...
    'VariableNames', {'Cluster1_PID','Cluster1_EF','Cluster2_PID','Cluster2_EF'});

%K means: Cluster 1 and Cluster 2 side by side

KM1_PID = KMClust1_PIDs;
KM1_EF  = EFKMnsClust_1;
KM2_PID = KMClust2_PIDs;
KM2_EF  = EFKMnsClust_2;

nKM = max(length(KM1_PID), length(KM2_PID));

KM1_PID_pad = [KM1_PID; nan(nKM-length(KM1_PID),1)];
KM1_EF_pad  = [KM1_EF;  nan(nKM-length(KM1_EF),1)];
KM2_PID_pad = [KM2_PID; nan(nKM-length(KM2_PID),1)];
KM2_EF_pad  = [KM2_EF;  nan(nKM-length(KM2_EF),1)];

KM_Table = table(KM1_PID_pad, KM1_EF_pad, KM2_PID_pad, KM2_EF_pad, ...
    'VariableNames', {'Cluster1_PID','Cluster1_EF','Cluster2_PID','Cluster2_EF'});

%Table 3: EF mean summary

SummaryTable = table( ...
    ["Hierarchical"; "Hierarchical"; "K means"; "K means"], ...
    [1; 2; 1; 2], ...
    [mean(EFHCClust_1); mean(EFHCClust_2); mean(EFKMnsClust_1); mean(EFKMnsClust_2)], ...
    'VariableNames', {'Method','Cluster','MeanEF'});

disp(SummaryTable)

%Save all three to Excel
outputfolder = pwd;
filename = fullfile(outputfolder,strcat('Clustering Results', file_name, sheet_name, '.xlsx'));
writetable(HC_Table,      filename, 'Sheet', 'Hierarchical')
writetable(KM_Table,      filename, 'Sheet', 'KMeans')
writetable(SummaryTable,  filename, 'Sheet', 'Summary')


% Box and Whisker Plots of Cluster Parameters

