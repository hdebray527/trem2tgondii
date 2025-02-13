data = readtable('GO bubble 021025.xlsx', 'ReadVariableNames', false); % Read without assuming headers
% Display first few rows to inspect
disp(data(1:5, :));
% Manually set correct column names
data.Properties.VariableNames = {'Pathway', 'GeneNumber', 'FDR', 'adjusted_P_value', 'logP'};
pathways = data.Pathway;
geneNumbers = data.GeneNumber;
fdrValues = data.FDR;
logPValues = data.logP;
genesList = data.Genes;

% Create a custom colormap: Red -> Pink -> Blue
num_colors = 256; % Number of color steps
red = [1, 0, 0];   % Red
pink = [1, 0.6, 0.6]; % Pink
blue = [0, 0, 1];   % Blue
% Blend the colors
colormap_custom = [linspace(red(1), pink(1), num_colors/2)', ...
linspace(red(2), pink(2), num_colors/2)', ...
linspace(red(3), pink(3), num_colors/2)'; ...
linspace(pink(1), blue(1), num_colors/2)', ...
linspace(pink(2), blue(2), num_colors/2)', ...
linspace(pink(3), blue(3), num_colors/2)'];
% Prepare the data
x = data.logP; % x-axis is now the logP values
y = 1:length(data.Pathway); % y-axis will be the indices of pathways
sizes = cellfun(@(x) numel(strsplit(x, ',')), data.Genes) * 27; % Bubble size: number of genes * 27
% Normalize FDR for color mapping
fdr = data.FDR; % FDR values for color
% Create bubble plot
figure;
scatter(x, y, sizes, fdr, 'filled'); % Color by FDR
colorbar; % Add color bar for FDR values
% Customize plot
set(gca, 'ytick', 1:length(data.Pathway), 'yticklabel', data.Pathway); % Set pathway names on y-axis
xlabel('logP');
ylabel('Pathway');
title('Bubble Plot of Pathways with logP and FDR');
% Apply the custom colormap
colormap(colormap_custom); % Use custom red-pink-blue colormap
% Adjust plot limits
xlim([min(x) * 0.9, max(x) * 1.1]);
ylim([0.5, length(data.Pathway) + 0.5]);
