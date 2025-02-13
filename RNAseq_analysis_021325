%% RNA sequencing analysis. Hannah, TREM2 KO vs WT
close all;
clear; clc;

% Load count table
count=readtable('counts_monocytes.txt');
sample_names={'KO1','KO2','KO3','KO4','WT1','WT2','WT3'};

% Gene names and values
gene_names=count.Geneid;
gene_length=count.Length/1000;
transcripts=count{:,7:end};

% Remove empty rows
remove_genes=find(sum(transcripts,2)==0);
transcripts(remove_genes,:)=[];
gene_names(remove_genes)=[];
gene_length(remove_genes)=[];

% Remove lowly expressed genes (less than 10 counts max)
min_counts=10;
remove_genes=find(max(transcripts,[],2)<min_counts);
transcripts(remove_genes,:)=[];
gene_names(remove_genes)=[];
gene_length(remove_genes)=[];

% Mitochondrial genes
mtgenes=contains(gene_names,'mt-');    

% Calculate TPM
RPK=transcripts./gene_length;
scaling_factor=sum(RPK)/1e6;
TPM=RPK./scaling_factor;

data=TPM;

%% PCA host genes (Supp. Fig. S2 C)
[coeff,score,latent,tsquared,explained] = pca(data');
f=figure; hold on;
for ii=1:4
    scatter(score(ii,1),score(ii,2),100,'filled','DisplayName',sample_names{ii},'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerFaceColor','w')
    text(score(ii,1)+500,score(ii,2),sample_names{ii})
end
for ii=5:7
    scatter(score(ii,1),score(ii,2),100,'filled','DisplayName',sample_names{ii},'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerFaceColor','k')
    text(score(ii,1)+500,score(ii,2),sample_names{ii})
end
xticks(''); yticks('');
xlabel ('PCA 1 (65%)'); ylabel ('PCA 2 (17%)')
set(gca,'FontSize',14)
f.Position=[680   558   560   420];
saveas(f,'PCA.jpg');

%% Differential Gene Expression (DGE)
T4 = array2table(upper(gene_names));
T5 = array2table(transcripts, 'VariableNames', sample_names);
T6 = [T4, T5];
writetable(T6, 'bulk_RNA.txt', 'Delimiter', '\t');

diffTable = rnaseqde(T6, 2:5, 6:8, IDColumns="Var1");
writetable(diffTable, 'KO_vs_WT_123124.xlsx');

% Volcano plot preparation (Fig. 4A)
pvals = [diffTable.AdjustedPValue];
log2FC = [diffTable.Log2FoldChange];

% Flip log2FC values to correct KO vs WT direction
log2FC = -log2FC;

% Define thresholds
log2FC_threshold = 1;
pval_threshold = 0.05;

% Identify differentially expressed genes
sig_up = pvals < pval_threshold & log2FC > log2FC_threshold;  % Upregulated
sig_down = pvals < pval_threshold & log2FC < -log2FC_threshold; % Downregulated

% Create figure for volcano plot
figure; hold on;

% Scatter plots
scatter(log2FC, -log10(pvals), 20, 'k', 'filled'); % Black dots, size 20
scatter(log2FC(sig_up), -log10(pvals(sig_up)), 40, 'r', 'filled'); % Red dots (Upregulated)
scatter(log2FC(sig_down), -log10(pvals(sig_down)), 40, 'b', 'filled'); % Blue dots (Downregulated)

% Add vertical dashed lines at Log2 Fold-Change = 1 and -1
xline(1, '--k', 'LineWidth', 1.5); % Black dashed line at Log2 Fold-Change = 1
xline(-1, '--k', 'LineWidth', 1.5); % Black dashed line at Log2 Fold-Change = -1

% Add horizontal dashed line at -log10 Adjusted p-value = 1.301
yline(1.30102999566, '--k', 'LineWidth', 1.5); % Black dashed line

% Annotate top significant genes
[~, topIdx] = sort(pvals); % Sort by p-value (ascending)
topGenes = topIdx(1:10); % Take top 10
x_shift = 0.1;  % Shift x position slightly
y_shift = 0.2;  % Shift y position slightly
text(log2FC(topGenes) + x_shift, -log10(pvals(topGenes)) + y_shift, gene_names(topGenes), ...
    'FontSize', 7, 'HorizontalAlignment', 'right');

% Labels and title
xlabel('Log2 Fold-Change (KO vs WT)');
ylabel('-log10 Adjusted p-value');
title('KO vs WT - Volcano Plot');

hold on; % Ensures subsequent plots can be added to the same figure

%% Table of TPM for significant hits only
% Ensure the correct number of columns for the data table
T10=array2table(gene_names(pvals<0.05 & (log2FC>1 | log2FC<-1)));
T11=array2table(data(pvals<0.05 & (log2FC>1 | log2FC<-1),:), 'VariableNames', sample_names);
T12=[T10, T11];
writetable(T12,'DGE_hits.txt','Delimiter','\t');

 %heat map (Fig. 4B)     
clustergram(data(pvals<0.05 & (log2FC>1 | log2FC<-1),:),'Standardize','row','Colormap',redbluecmap);
