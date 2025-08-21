%% Fig. 2a, Raster & FR
clear -path; 
dataName = 'fig2a_S13_2_HPC_chan_13_cluster_1'; YLim = [-2 10]; YTickRange = [-2 0 10]; YTickLabel = [-2 0 10];
load([path '\fig2\' dataName '.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
fontNumbleLable = 5; fontNumbleTick = 5;
LineWidth = 0.25; 
XLim = [-1 3]; XTickRange = [-1:3]; XTickLabel = [-1:3];
close all;
PlotExampleNeuronalActivity;

outputFigure(fullfile([path dataName(1:4)]), dataName);

%% Fig. 2a, spike waveform
clear -path; 
fontNumbleLable = 6; fontNumbleTick = 6;
LineWidth = 0.25;

close all;
load(['data_spikes.mat']);
load(['times_data.mat']);
idx = find(cluster_class(:, 1) == clusters_check(i));
f = figure('units', 'centimeters', 'position', [5, 10, 1.5, 2], 'Color', [1, 1, 1]);
ax5 =   axes('units', 'centimeters', 'position', [0.5, .5, 1, 1]); hold on;
outPlot = LK_SpikeDensityPlot(spikes(idx, :), par.sr);
colormap(ax5, 'parula');
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
    'xlim', round([min(outPlot.spikeTime), max(outPlot.spikeTime)]), 'xtick', round([min(outPlot.spikeTime), max(outPlot.spikeTime)]), ...
    'ylim', [outPlot.yMin, outPlot.yMax], 'ytick', [outPlot.yMin, outPlot.yMax], 'ticklength', [0, 0]);
xl = xlabel('ms', 'position', [1, outPlot.yMin - (outPlot.yMax - outPlot.yMin) * 0.1, -1]);
yl = ylabel('\muV', 'position', [-0.35, (outPlot.yMax + outPlot.yMin)/2, -1]);
% indicate number of spikes used for the analysis
tx = text(1, outPlot.yMax + (outPlot.yMax - outPlot.yMin) * 0.1, ['n=', num2str(length(idx))], 'HorizontalAlignment', 'center', 'FontSize', fontNumbleTick);
outputFigure(fullfile([path 'fig2']), ['SpikeWaveform']);


%% Fig. 2b/c, population 
clear -path; 
dataName = 'fig2b_AMY';
% dataName = 'fig2b_HPC';
load([path '\fig2\' dataName '.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
fontNumbleLable = 6; fontNumbleTick = 6;
LineWidth = 0.25; 
XLim = [-1 3]; XTickRange = [-1:3]; XTickLabel = [-1:3];
YLim = [-0.3 1; -0.8 1]; YTickRange = [-0.3 0 1; -0.8 0 1]; YTickLabel = [-0.3 0 1; -0.8 0 1];
PlotPopulationActivity;


%% Fig. 2e/f, GLM 
clear -path; 
% dataName = 'fig2b_AMY';
dataName = 'fig2b_HPC';
load([path '\fig2\' dataName '.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
fontNumbleLable = 6; fontNumbleTick = 6;
LineWidth = 0.25; 
XLim = [-1 3]; XTickRange = [-1:3]; XTickLabel = [-1:3];
YLim = [-0.3 1; -0.8 1]; YTickRange = [-0.3 0 1; -0.8 0 1]; YTickLabel = [-0.3 0 1; -0.8 0 1];
PlotPopulationGLM


