function Cluster(file_name, sheet_name, Norm_Flag, EF_Status, Num_Clusts)

close all

% The Cluster function clusters patients based on their parameters via 
% hierarchical and k means clustering methods

% Input the name of your file and specify the sheet you want data read from
% Input as a string ie, 'Patient Optimized Parameters','Sheet1'

% Norm_Flag determines what type of data normalization is conducted 
% 0 for no normalization, 1 for normalization by the mean, and 2 for normalization by the standard deviation

%EF_Status determines whether or not you wish to include Ejection Fraction
%(EF) as a metric by which the algorithms will cluster patients

% 0 for no and 1 for yes

% Set Num_Clusts to the number of clusters you want, 2 or 3

% Goes back into directory to find spreadsheet
cd ..

% This function applies both hierarchical and k means clustering to your
% input data set

rng(1);

% Input your data set as a string, EX: Cluster('Optimized parameters')

data_table = readtable(file_name,'Sheet',sheet_name,'VariableNamingRule', 'preserve', 'NumHeaderLines', 1, 'ReadVariableNames', true);      

% Extracts list of patient IDs from 1st row excel header 

PatIDs = fillmissing(string(readcell(file_name, 'Range', '1:1')), 'previous');  % reads first excel row into strings filling empty cells with previous value
patient_strings = unique(PatIDs(contains(PatIDs, 'Patient')), 'stable'); % Filters to extract PIDs and preservers order
PIDs = str2double(regexp(patient_strings, '\d+', 'match', 'once')); % finds the digits and turns the string into a PID vector
PIDs = PIDs';
PIDs_Index = containers.Map(PIDs, 1:numel(PIDs));

 
% Getting corresponding PID Ejection Fractions (EF) from patient data 
cd ../'CV Model'/Core/
patientdata_table = readtable('AllPatsMRI_2_REU.xlsx','Sheet','Sheet1','VariableNamingRule','preserve');
[~, PID_indices] = ismember(PIDs, patientdata_table.Pat_Num);
EF = patientdata_table.EF(PID_indices);

% Calculating mPAP
PpaS = patientdata_table.PpaS(PID_indices);
PpaD = patientdata_table.PpaD(PID_indices);
PsaS = patientdata_table.PsaS(PID_indices);
PsaD = patientdata_table.PsaD(PID_indices);

mPAP = (PpaS + (2*PpaD))/3;
MAP = (PsaS + (2*PsaD))/3;

PH_patients = mPAP > 20;

cd ..
cd ../'CV Model'/Clustering/

% Preparing data for normalization
if (EF_Status == 0)
mean_column_flags = contains(data_table.Properties.VariableNames, 'Mean');
mean_matrix_raw = table2array(data_table(:,mean_column_flags));
A_Optp = mean_matrix_raw';
[NumPats_Optp, Num_Optp] = size(A_Optp);
elseif (EF_Status == 1)
mean_column_flags = contains(data_table.Properties.VariableNames, 'Mean');
mean_matrix_raw = table2array(data_table(:,mean_column_flags));
A_Optp = mean_matrix_raw';
A_Optp = [A_Optp, EF];
[NumPats_Optp, Num_Optp] = size(A_Optp);
end 

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

% 1. Calculate cumulative variance percentage (%)
var_explained = (sigma_Optp.^2) / sum(sigma_Optp.^2) * 100;
cum_var_explained = cumsum(var_explained);

% 2. Create Plot
figure('Name', 'Cumulative Variance Explained', 'Position', [100 100 700 450]);

plot(1:length(cum_var_explained), cum_var_explained, '-o', ...
    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.8], 'Color', [0.2 0.5 0.8]);

% Formatting
title('Cumulative Variance Explained by Principal Components', 'FontSize', 14);
xlabel('Principal Component Number', 'FontSize', 12);
ylabel('Cumulative Variance Explained (%)', 'FontSize', 12);
xticks(1:length(cum_var_explained));
ylim([0 105]);
grid on;


%Perform hierarchical clustering and make a dendrogram

% Creating hierarchical tree where we use the Ward metric when joining clusters which 
%  considers the increase in the within cluster sum of squares. This is just the squared 
%  distance between each element in the cluster and it's centroid
NumHCClust_Optp = Num_Clusts;
HCLink_Optp = linkage(ANorm_Optp,'ward');                   % Making the tree
HCClust_Optp = cluster(HCLink_Optp, ...                     % Cluster number to display
'Maxclust',NumHCClust_Optp); 
if (Num_Clusts == 2)
HCCutOff_Optp = median([HCLink_Optp(end-1,3) HCLink_Optp(end,3)]);
elseif (Num_Clusts == 3)
    HCCutOff_Optp = median([HCLink_Optp(end-1,3) HCLink_Optp(end-2,3)]);
end 

%Performing kmeans
% Running KMeans clustering on normalized optimized parameters

KMnsClust_Num = Num_Clusts;

[KMnsClust_Optp,KMnsCntrd_Optp, ...                 % KMeans clust numbers
KMnsSumDst_Optp] = ...                          %  centroids, sum of dist
kmeans(ANorm_Optp,KMnsClust_Num, ...            %  Norm data, # of clusts 
'Distance','cityblock', ...                     %  Use L1 distance
'Display','off', ...                          %  Display results
'Replicates',20);                               %  Run 20 times  

% Cluster based on Ejection Fraction (EF)

EFclass = zeros(size(EF));
EFclass(EF < 0.5) = 1;
EFclass(EF >= 0.5) = 2;


% Separating HFpEF pateints from ANorm_Opt
HFpEF = ANorm_Optp(EFclass == 2,:);
HFrEF = ANorm_Optp(EFclass == 1,:);

% K means and Hierarchical Cluster Alignment - Claude 

C = confusionmat(HCClust_Optp, KMnsClust_Optp);
perms_list = perms(1:Num_Clusts);
best_agree = -1;
best_perm = 1:Num_Clusts;
for p = 1:size(perms_list,1)
    agree_count = sum(diag(C(:, perms_list(p,:))));
    if agree_count > best_agree
        best_agree = agree_count;
        best_perm = perms_list(p,:);
    end
end
map = zeros(Num_Clusts,1);
map(best_perm) = 1:Num_Clusts;
KMnsClust_Optp = map(KMnsClust_Optp);

% Canonicalize cluster labels using mean EF, so labels/colors are 
% comparable across all runs and configurations (Num_Clusts, EF_Status)
clust_meanEF = zeros(Num_Clusts,1);
for k = 1:Num_Clusts
    clust_meanEF(k) = mean(EF(HCClust_Optp == k));
end
[~, order] = sort(clust_meanEF, 'ascend');   % lowest EF -> label 1
canonical_map = zeros(Num_Clusts,1);
canonical_map(order) = 1:Num_Clusts;

HCClust_Optp   = canonical_map(HCClust_Optp);
KMnsClust_Optp = canonical_map(KMnsClust_Optp);

% Getting colors 
cmap = parula(5);       % sample 5 points across parula for good spread
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
if (Num_Clusts == 2)
tab2 = uitab(tg, 'Title', 'PCA - Hierarchical');
ax2 = axes('Parent', tab2);
gscatter(ax2, AV_Optp(:,1), AV_Optp(:,2), HCClust_Optp, [c1; c3], '.', 35)
set(ax2, 'XTick', [], 'YTick', [], 'FontSize', 16)
xlabel(ax2, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp), 'FontSize', 18)
ylabel(ax2, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp), 'FontSize', 18)
title(ax2, 'PCA of Optimized Parameters, Clustered via Hierarchical', 'FontSize', 18)
legend(ax2, 'Cluster 1', 'Cluster 2', 'FontSize', 14)
elseif (Num_Clusts == 3)
    tab2 = uitab(tg, 'Title', 'PCA - Hierarchical');
    ax2 = axes('Parent', tab2);
    gscatter(ax2, AV_Optp(:,1), AV_Optp(:,2), HCClust_Optp, [c1; c3; c4], '.', 35)
    set(ax2, 'XTick', [], 'YTick', [], 'FontSize', 16)
    xlabel(ax2, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp), 'FontSize', 18)
    ylabel(ax2, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp), 'FontSize', 18)
    title(ax2, 'PCA of Optimized Parameters, Clustered via Hierarchical', 'FontSize', 18)
    legend(ax2, 'Cluster 1', 'Cluster 2','Cluster 3','FontSize', 14)
end 


color = [c1; c1; c3; c3];
symbol = '.x.x';

%Tab 3: KMeans PCA scatter
if (Num_Clusts == 2)
tab3 = uitab(tg, 'Title', 'PCA - KMeans');
ax3 = axes('Parent', tab3);
h3 = gscatter(ax3, AV_Optp(:,1),AV_Optp(:,2), {KMnsClust_Optp, PH_patients},color,symbol,35);
h3(2).MarkerSize = 15;   % Cluster 1, PH — the 'x' markers
h3(4).MarkerSize = 15;   % Cluster 2, PH — the 'x' markers
set(ax3, 'XTick', [], 'YTick', [],'FontSize',16)
xlabel(ax3, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp),'FontSize', 22)
ylabel(ax3, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp),'FontSize', 22)
title(ax3, 'PCA of Optimized Parameters, Clustered via k Means', 'FontSize', 18)
legend(ax3, {'Cluster 1, Non-PH','Cluster 1, PH','Cluster 2, Non-PH','Cluster 2, PH'}, 'FontSize', 14)
elseif (Num_Clusts == 3)
    tab3 = uitab(tg, 'Title', 'PCA - KMeans');
    ax3 = axes('Parent', tab3);
    gscatter(ax3, AV_Optp(:,1),AV_Optp(:,2),KMnsClust_Optp,[c1; c3; c4],'.',35)
    set(ax3, 'XTick', [], 'YTick', [],'FontSize',16)
    xlabel(ax3, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp),'FontSize', 18)
    ylabel(ax3, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp),'FontSize', 18)
    title(ax3, 'PCA of Optimized Parameters, Clustered via k Means', 'FontSize', 18)
    legend(ax3, 'Cluster 1', 'Cluster 2','Cluster 3','FontSize',14)
end 

%Tab 3b: KMeans PCA scatter, marker shape = HFrEF / HFpEF
tab3b = uitab(tg, 'Title', 'PCA - KMeans (EF Class)');
ax3b = axes('Parent', tab3b);

if (Num_Clusts == 2)
    colorEF  = [c1; c1; c3; c3];
    symbolEF = '.x.x';
    h3b = gscatter(ax3b, AV_Optp(:,1), AV_Optp(:,2), {KMnsClust_Optp, EFclass}, colorEF, symbolEF, 35);
    h3b(2).MarkerSize = 15;   % Cluster 1, HFpEF 'x' markers
    h3b(4).MarkerSize = 15;   % Cluster 2, HFpEF 'x' markers
    legend(ax3b, {'Cluster 1, HFrEF','Cluster 1, HFpEF','Cluster 2, HFrEF','Cluster 2, HFpEF'}, 'FontSize', 14)
elseif (Num_Clusts == 3)
    colorEF  = [c1; c1; c3; c3; c4; c4];
    symbolEF = '.x.x.x';
    h3b = gscatter(ax3b, AV_Optp(:,1), AV_Optp(:,2), {KMnsClust_Optp, EFclass}, colorEF, symbolEF, 35);
    h3b(2).MarkerSize = 15;
    h3b(4).MarkerSize = 15;
    h3b(6).MarkerSize = 15;
    legend(ax3b, {'Cluster 1, HFrEF','Cluster 1, HFpEF','Cluster 2, HFrEF','Cluster 2, HFpEF','Cluster 3, HFrEF','Cluster 3, HFpEF'}, 'FontSize', 14)
end

set(ax3b, 'XTick', [], 'YTick', [], 'FontSize', 16)
xlabel(ax3b, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp), 'FontSize', 22)
ylabel(ax3b, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp), 'FontSize', 22)
title(ax3b, 'PCA of Optimized Parameters, Clustered via k Means', 'FontSize', 18)

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

if (Num_Clusts == 2)
    ThreeGroup = zeros(size(HCClust_Optp));
    ThreeGroup(HCClust_Optp == 1 & KMnsClust_Optp == 1) = 1;
    ThreeGroup(HCClust_Optp == 2 & KMnsClust_Optp == 2) = 2;
    ThreeGroup(HCClust_Optp ~= KMnsClust_Optp) = 3;

    gscatter(ax7, AV_Optp(:,1), AV_Optp(:,2), ThreeGroup, [c1; c3; c5], '.', 35)
    legend(ax7, 'Agreed Cluster 1','Agreed Cluster 2','Disagree','FontSize', 14)

elseif (Num_Clusts == 3)
    ThreeGroup = zeros(size(HCClust_Optp));
    ThreeGroup(HCClust_Optp == 1 & KMnsClust_Optp == 1) = 1;
    ThreeGroup(HCClust_Optp == 2 & KMnsClust_Optp == 2) = 2;
    ThreeGroup(HCClust_Optp == 3 & KMnsClust_Optp == 3) = 3;
    ThreeGroup(HCClust_Optp ~= KMnsClust_Optp) = 4;

    gscatter(ax7, AV_Optp(:,1), AV_Optp(:,2), ThreeGroup, [c1; c3; c4; c5], '.', 35)
    legend(ax7, 'Agreed Cluster 1','Agreed Cluster 2','Agreed Cluster 3','Disagree','FontSize', 14)
end

 num_disagree = sum(ThreeGroup == 3);
    pct_disagree = (num_disagree / numel(ThreeGroup)) * 100;
    
    fprintf('Total disagreeing points: %d out of %d (%.2f%%)\n', ...
        num_disagree, numel(ThreeGroup), pct_disagree);

set(ax7, 'XTick', [], 'YTick', [],'FontSize',16)
xlabel(ax7, sprintf('PC1 (%.1f%% variance)', 100*sigma_Optp(1)^2/rho_Optp),'FontSize', 18)
ylabel(ax7, sprintf('PC2 (%.1f%% variance)', 100*sigma_Optp(2)^2/rho_Optp), 'FontSize', 18)
title(ax7, 'Agreement Between Hierarchical and K Means Clustering','FontSize', 18)

% Tab 8: PH Confusion Matrix
tab8 = uitab(tg, 'Title', 'Pulmonary Hypertension Confusion Matrix');

clustLabels = categorical(KMnsClust_Optp, [1 2], {'Cluster 1','Cluster 2'});
phLabels = categorical(PH_patients, [0 1], {'Non-PH','PH'});

cm8 = confusionchart(tab8, clustLabels, phLabels);
title(cm8, 'Cluster vs PH Status')

% Tab 9: HFpEF/HFrEF Confusion Matrix
tab9 = uitab(tg, 'Title', 'HF Status Confusion Matrix');
clustLabels = categorical(KMnsClust_Optp, [1 2], {'Cluster 1','Cluster 2'});
phLabels = categorical(EFclass, [1 2], {'HFrEF','HFpEF'});

cm9 = confusionchart(tab9, clustLabels, phLabels);
title(cm9, 'Cluster vs HF Status')

% Tab 10: EF vs mPAP
tab10 = uitab(tg, 'Title', 'EF vs mPAP');
ax10 = axes('Parent', tab10);

scatter(ax10, EF, mPAP, 60, 'filled', 'MarkerFaceColor', c1);
hold(ax10, 'on')
yline(ax10, 20, '--k', 'mPAP = 20 mmHg', 'LabelHorizontalAlignment', 'left');   % PH threshold used above
xline(ax10, 0.5, '--k', 'EF = 0.5', 'LabelVerticalAlignment', 'bottom');        % HFrEF / HFpEF cutoff

set(ax10, 'FontSize', 16)
xlabel(ax10, 'Ejection Fraction (EF)', 'FontSize', 18)
ylabel(ax10, 'Mean Pulmonary Artery Pressure (mmHg)', 'FontSize', 18)
title(ax10, 'EF vs mPAP', 'FontSize', 18)
grid(ax10, 'on')
box(ax10, 'on')


% % Clustering Summary Table
% 
% Getting PIDs of each cluster
HCClust1_PIDs = PIDs(HCClust_Optp == 1);
HCClust2_PIDs = PIDs(HCClust_Optp == 2);
KMClust1_PIDs = PIDs(KMnsClust_Optp == 1);
KMClust2_PIDs = PIDs(KMnsClust_Optp == 2);
KMClust3_PIDs = PIDs(KMnsClust_Optp == 3);
EFHCClust_1 = EF(HCClust_Optp == 1);
EFHCClust_2 = EF(HCClust_Optp == 2);
EFKMnsClust_1 = EF(KMnsClust_Optp == 1);
EFKMnsClust_2 = EF(KMnsClust_Optp == 2);
MeanEF_HC1 = mean(EFHCClust_1);
MeanEF_HC2 = mean(EFHCClust_2);
MeanEF_KM1 = mean(EFKMnsClust_1);
MeanEF_KM2 = mean(EFKMnsClust_2);
% 
% % Creating Table of PIDs, Cluster, and EF - Claude
% 
% %Hierarchical: Cluster 1 and Cluster 2 side by side
% 
% HC1_PID = HCClust1_PIDs;
% HC1_EF  = EFHCClust_1;
% HC2_PID = HCClust2_PIDs;
% HC2_EF  = EFHCClust_2;
% 
% nHC = max(length(HC1_PID), length(HC2_PID));
% 
% HC1_PID_pad = [HC1_PID; nan(nHC-length(HC1_PID),1)];
% HC1_EF_pad  = [HC1_EF;  nan(nHC-length(HC1_EF),1)];
% HC2_PID_pad = [HC2_PID; nan(nHC-length(HC2_PID),1)];
% HC2_EF_pad  = [HC2_EF;  nan(nHC-length(HC2_EF),1)];
% 
% HC_Table = table(HC1_PID_pad, HC1_EF_pad, HC2_PID_pad, HC2_EF_pad, ...
%     'VariableNames', {'Cluster1_PID','Cluster1_EF','Cluster2_PID','Cluster2_EF'});
% 
% %K means: Cluster 1 and Cluster 2 side by side
% 
% KM1_PID = KMClust1_PIDs;
% KM1_EF  = EFKMnsClust_1;
% KM2_PID = KMClust2_PIDs;
% KM2_EF  = EFKMnsClust_2;
% 
% nKM = max(length(KM1_PID), length(KM2_PID));
% 
% KM1_PID_pad = [KM1_PID; nan(nKM-length(KM1_PID),1)];
% KM1_EF_pad  = [KM1_EF;  nan(nKM-length(KM1_EF),1)];
% KM2_PID_pad = [KM2_PID; nan(nKM-length(KM2_PID),1)];
% KM2_EF_pad  = [KM2_EF;  nan(nKM-length(KM2_EF),1)];
% 
% KM_Table = table(KM1_PID_pad, KM1_EF_pad, KM2_PID_pad, KM2_EF_pad, ...
%     'VariableNames', {'Cluster1_PID','Cluster1_EF','Cluster2_PID','Cluster2_EF'});
% 
% %Table 3: EF mean summary
% 
% SummaryTable = table( ...
%     ["Hierarchical"; "Hierarchical"; "K means"; "K means"], ...
%     [1; 2; 1; 2], ...
%     [mean(EFHCClust_1); mean(EFHCClust_2); mean(EFKMnsClust_1); mean(EFKMnsClust_2)], ...
%     'VariableNames', {'Method','Cluster','MeanEF'});
% 
% 
% %Save all three to Excel
% outputfolder = pwd;
% filename = fullfile(outputfolder,'Clustering Results.xlsx');
% writetable(HC_Table,      filename, 'Sheet', 'Hierarchical')
% writetable(KM_Table,      filename, 'Sheet', 'KMeans')
% writetable(SummaryTable,  filename, 'Sheet', 'Summary')

%Tab 8: Box and Whisker Plots by Cluster - Hierarchical

cd ..
paramNames = readcell(file_name, 'Sheet', sheet_name, 'Range', 'B3:B11');
if (EF_Status == 1)
paramNames = [paramNames; 'EF'];
end

ElvM = A_Optp(:,6);
Elvm = A_Optp(:,7);
ErvM = A_Optp(:,8);
Ervm = A_Optp(:,9);

Elv = ElvM - Elvm;
Erv = ErvM - Ervm;

params = [A_Optp, EF, Elv, Erv, mPAP, MAP];
paramNames_EF = [paramNames; 'EF'; 'ElvM - Elvm'; 'ErvM - Ervm'; 'mPAP';'MAP'];
[~, Num_Optp_EF] = size(params);

figure;
t = tiledlayout('flow');
title(t, 'Parameter Values Across Patients - Hierarchical');
if (Num_Clusts == 2)
for c = 1:Num_Optp
    nexttile;
    HCClust1_params = A_Optp(HCClust_Optp == 1, c);
    HCClust2_params = A_Optp(HCClust_Optp == 2, c);
    combinedData = [HCClust1_params; HCClust2_params];
    groupLabels = [repmat({'Cluster 1'}, numel(HCClust1_params), 1); repmat({'Cluster 2'}, numel(HCClust2_params), 1)];
    b = boxplot(combinedData, groupLabels, 'Widths', 0.6);
    set(b, 'LineWidth', 2);
    ylabel(paramNames(c));
    grid on;
    

end
elseif (Num_Clusts == 3)
for c = 1:Num_Optp
    nexttile;
    HCClust1_params = A_Optp(HCClust_Optp == 1, c);
    HCClust2_params = A_Optp(HCClust_Optp == 2, c);
    HCClust3_params = A_Optp(HCClust_Optp == 3, c);
    combinedData = [HCClust1_params; HCClust2_params; HCClust3_params];
    groupLabels = [repmat({'Cluster 1'}, numel(HCClust1_params), 1); repmat({'Cluster 2'}, numel(HCClust2_params), 1);repmat({'Cluster 3'}, numel(HCClust3_params), 1)];
    b = boxplot(combinedData, groupLabels, 'Widths', 0.6);
    set(b, 'LineWidth', 2);
    ylabel(paramNames(c));
    grid on;
end
end 

% % Check for statistical differences
% % First perform t-test
% disp('Results from t-test for each parameter:')
% for i = 1:numel(paramNames)
%     KMnsClust1_params = A_Optp(KMnsClust_Optp == 1, i);
%     KMnsClust2_params = A_Optp(KMnsClust_Optp == 2, i);
%     [h, p, ci, stats] = ttest2(KMnsClust1_params(:,:), KMnsClust2_params(:,:));
%     fprintf('The p-value for %s is %.4f. The h-value is %d.\n', paramNames{i}, p, h);
% end

figure;
t = tiledlayout('flow');
title(t, 'Parameter Values Across Patients - K Means');
if (Num_Clusts == 2)
for c = 1:Num_Optp_EF
    nexttile;
    KMnsClust1_params = params(KMnsClust_Optp == 1, c);
    KMnsClust2_params = params(KMnsClust_Optp == 2, c);
    combinedData = [KMnsClust1_params; KMnsClust2_params];
    groupLabels = [repmat({'Cluster 1'}, numel(KMnsClust1_params), 1); repmat({'Cluster 2'}, numel(KMnsClust2_params), 1)];
    b = boxplot(combinedData, groupLabels, 'Widths', 0.6);
    set(b, 'LineWidth', 2);
    ylabel(paramNames_EF(c));
    grid on;
    [p] = ttest2(KMnsClust1_params(:,:), KMnsClust2_params(:,:));
    % p_str = sprintf('p = %.2e', p);
    % text(0.5, 0.90, p_str, 'Units', 'normalized', ...
    % 'HorizontalAlignment', 'center', ...
    % 'FontSize', 10, ...
    % 'FontWeight', 'bold');
end
elseif (Num_Clusts == 3)
for c = 1:Num_Optp
    nexttile;
    KMnsClust1_params = params(KMnsClust_Optp == 1, c);
    KMnsClust2_params = params(KMnsClust_Optp == 2, c);
    KMnsClust3_params = params(KMnsClust_Optp == 3, c);
    combinedData = [KMnsClust1_params; KMnsClust2_params; KMnsClust3_params];
    groupLabels = [repmat({'Cluster 1'}, numel(KMnsClust1_params), 1); repmat({'Cluster 2'}, numel(KMnsClust2_params), 1);repmat({'Cluster 3'}, numel(KMnsClust3_params), 1)];
    b = boxplot(combinedData, groupLabels, 'Widths', 0.6);
    set(b, 'LineWidth', 2);
    ylabel(paramNames(c));
    grid on;
    % [p, h] = ranksum(KMnsClust1_params(:,:), KMnsClust2_params(:,:), KMnsClust3_params(:,:));
    % p_str = sprintf('p = %.2e', p);
    % text(0.5, 0.90, p_str, 'Units', 'normalized', ...
    % 'HorizontalAlignment', 'center', ...
    % 'FontSize', 10, ...
    % 'FontWeight', 'bold');
end
end 

% Parameter Distributions

figure;
t = tiledlayout('flow');
title(t, 'Parameter Distributions Across Patients - K Means');
if (Num_Clusts == 2)
for c = 1:Num_Optp_EF
    nexttile;
    KMnsClust1_params = params(KMnsClust_Optp == 1, c);
    KMnsClust2_params = params(KMnsClust_Optp == 2, c);
    [y1,x1] = kde(KMnsClust1_params);
    [y2,x2] = kde(KMnsClust2_params);

    [~, bin_edges] = histcounts(params(:, c), 25);

    y1_norm = y1*mean(diff(bin_edges));
    y2_norm = y2*mean(diff(bin_edges));

    [p, h] = ranksum(KMnsClust1_params, KMnsClust2_params);

    % histogram(KMnsClust1_params, bin_edges, 'Normalization','probability', 'FaceColor', c1, 'FaceAlpha', 0.5, 'EdgeColor','none')
    hold on
    hold on
    plot(x1, y1_norm, '-', 'Color', c1, 'LineWidth', 2)
    % histogram(KMnsClust2_params, bin_edges, 'Normalization','probability', 'FaceColor', c3, 'FaceAlpha', 0.5, 'EdgeColor','none')
    hold on 
    plot(x2, y2_norm, '-', 'Color', c3, 'LineWidth', 2)
    ylabel('Normalized Probability')
    xlabel(paramNames_EF(c),'FontSize',16);
    legend({'Cluster 1','Cluster 2'}, 'FontSize', 8)
    grid on;

    % p_str = sprintf('p = %.2e', p);
    % 
    % text(0.5, 0.90, p_str, 'Units', 'normalized', ...
    % 'HorizontalAlignment', 'center', ...
    % 'FontSize', 10, ...
    % 'FontWeight', 'bold');
    
end
elseif (Num_Clusts == 3)
    for c = 1:Num_Optp
    nexttile;
    KMnsClust1_params = params(KMnsClust_Optp == 1, c);
    KMnsClust2_params = params(KMnsClust_Optp == 2, c);
    KMnsClust3_params = params(KMnsClust_Optp == 3, c);
    [y1,x1] = kde(KMnsClust1_params);
    [y2,x2] = kde(KMnsClust2_params);
    [y3,x3] = kde(KMnsClust3_params);

    [~, bin_edges] = histcounts(params(:, c), 25);
    
    y1_norm = y1*mean(diff(bin_edges));
    y2_norm = y2*mean(diff(bin_edges));
    y3_norm = y3*mean(diff(bin_edges));
    
    % [p h] = ranksum(KMnsClust1_params, KMnsClust2_params);
    
    histogram(KMnsClust1_params, bin_edges, 'Normalization','probability', 'FaceColor', c1, 'FaceAlpha', 0.5, 'EdgeColor','none')
    hold on
    hold on
    plot(x1, y1_norm, '-', 'Color', c1, 'LineWidth', 2, 'HandleVisibility','off')
    histogram(KMnsClust2_params, bin_edges, 'Normalization','probability', 'FaceColor', c3, 'FaceAlpha', 0.5, 'EdgeColor','none')
    hold on 
    plot(x2, y2_norm, '-', 'Color', c2, 'LineWidth', 2, 'HandleVisibility','off')
    hold on
    histogram(KMnsClust3_params, bin_edges, 'Normalization', 'probability','FaceColor',c4,'FaceAlpha',0.5,'EdgeColor','none')
    hold on 
    plot(x3, y3_norm, '-', 'Color', c4*0.7,'LineWidth', 2, 'HandleVisibility', 'off')
    ylabel('Normalized Probability')
    xlabel(paramNames(c));
    legend({'Cluster 1','Cluster 2'}, 'FontSize', 8)
    grid on;
    
    % p_str = sprintf('p = %.2e', p);
    % 
    % text(0.5, 0.90, p_str, 'Units', 'normalized', ...
    % 'HorizontalAlignment', 'center', ...
    % 'FontSize', 10, ...
    % 'FontWeight', 'bold');
    end
end 

% Table representation of clutered parameters
outputFile = 'Clustered Parameters.xlsx';

for c = 1:Num_Optp_EF
    KMnsClust1_params = params(KMnsClust_Optp == 1, c);
    KMnsClust2_params = params(KMnsClust_Optp == 2, c);

    n1 = numel(KMnsClust1_params);
    n2 = numel(KMnsClust2_params);
    maxN = max(n1, n2);

    PID1 = NaN(maxN, 1);
    PID2 = NaN(maxN, 1);
    Val1 = NaN(maxN, 1);
    Val2 = NaN(maxN, 1);

    PID1(1:n1) = KMClust1_PIDs;
    PID2(1:n2) = KMClust2_PIDs;
    Val1(1:n1) = KMnsClust1_params;
    Val2(1:n2) = KMnsClust2_params;

    T = table(PID1, Val1, PID2, Val2, ...
        'VariableNames', {'PID_Cluster1','Value_Cluster1','PID_Cluster2','Value_Cluster2'});

    % Use plain char for the parameter name, not a cell/string extraction
     paramName = char(paramNames_EF(c));
    sheetName = sprintf('%d_%s', c, paramName(1:min(29,end)));  % e.g. "3_Elvm"
    writetable(T, outputFile, 'Sheet', sheetName);
end

cd Clustering/



