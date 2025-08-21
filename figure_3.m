%% Fig. 3a/b, time-frequency decomposition
addpath 'E:\Documents\Matlab Support\fieldtrip-20240603';
ft_defaults
clear -path; 
dataName = 'fig3a';
load([path '\fig3\' dataName '.mat' ])
fontNumbleLable = 6; fontNumbleTick = 6;
LineWidth = 0.25;
f_min_1 = 1; f_max_1 = 30; YTickRange = [1 10 20 30]; XTickRange = -1:3; XTickLabel = [-1:3]; freqType = 'low';
for i_reg = 1:n_reg
    %------ pow of all valence conditions
    close all;
    f = figure('units', 'centimeters', 'position', [5 5 4 2.5], 'Color', [1, 1, 1]);
    pow_all = [];
    for i = 1:length(trials_valence)
        pow_all = cat(1, pow_all, pow_valence{i_reg, i});
    end

    freq_tmp.dimord = dimord;
    freq_tmp.label = {'avg'};
    freq_tmp.powspctrm = nanmean(pow_all, dim);

    cfg = [];
    cfg.toi = 'all';
    cfg.xlim = [bl_window(1) ts_window(2)];
    cfg.ylim = [f_min_1 f_max_1];
    cfg.maskstyle = 'saturation';
    cfg.parameter = 'powspctrm';
    cfg.figure = 'gca';
    cfg.title = ' ';
    ft_singleplotTFR(cfg, freq_tmp);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    clim([-0.2 0.2]);

    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
            'YTick', YTickRange, 'YTickLabel', YTickRange,'XTick', XTickRange, 'XTickLabel', XTickRange, ...
            'tickdir', 'out', 'ticklength', [0.02, 0.02]); box off;
    c = colorbar; c.Ticks = [-0.2 0 0.2];
    
    %------------ permutation test
    freq_tmp.powspctrm = pow_all;
    cfg = [];
    cfg.latency = bl_window; % Baseline interval
    freq_baseline = ft_selectdata(cfg, freq_tmp);
    pow_tmp1 = freq_baseline.powspctrm;
    cfg.latency = [bl_window(1) ts_window(2)]; % Task interval
    freq_epoch = ft_selectdata(cfg, freq_tmp);
    pow_tmp2 = freq_epoch.powspctrm;
    pow_tmp1 = repmat(mean(pow_tmp1, 4), [1, 1, 1, size(pow_tmp2, 4)]);
    d1 = [];
    d1.time = freq_epoch.time;
    d1.freq = freq_epoch.freq;
    d1.label = {'avg'};
    d1.dimord = freq_epoch.dimord;
    d2 = d1;
    d1.powspctrm = pow_tmp1;
    d2.powspctrm = pow_tmp2;
    if strcmp(shuffle_dim, 'channel')
        [a, b, c, d] = size(pow_tmp1);
        d1.powspctrm = reshape(pow_tmp1, b, a, c, d);
        d2.powspctrm = reshape(pow_tmp2, b, a, c, d);
    end
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    cfg.latency             = [bl_window(1) ts_window(2)];
    cfg.frequency           = [f_min_1 f_max_1];
    cfg.statistic           = 'depsamplesT';
    cfg.alpha               = 0.05;
    cfg.clusteralpha        = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.tail                = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
    cfg.clustertail         = cfg.tail;
    cfg.clusterstatistic    = 'wcm';
    cfg.numrandomization    = 1000;
    cfg.neighbours          = [];
    n_min = size(pow_tmp1, dim);
    cfg.design(1, 1:2*n_min) = [ones(1, n_min) 2*ones(1, n_min)];
    cfg.design(2, 1:2*n_min) = [1:n_min 1:n_min];
    cfg.ivar = 1;
    cfg.uvar = 2;
    stat    = ft_freqstatistics(cfg, d1, d2);

    sig     = squeeze(stat.mask);
    idxF = find(stat.freq >= f_min_1 & stat.freq <= f_max_1)
    plot_contour(stat.time, stat.freq,double(sig), 'k-', LineWidth);
    

    outputFigure(fullfile([path '\fig3\']), [dataName '_' regions{i_reg} '_' freqType '_1']);
    
end

%-- Compare pow by averaging the frequency band then apply permutation test
addpath 'E:\Documents\Matlab Support\fieldtrip-20240603';
ft_defaults
addpath('E:\Projects\IAPS\Code\auxiliary');
clear -path;
dataName = 'fig3a';
load([path '\fig3\' dataName '.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
fontNumbleLable = 6; fontNumbleTick = 6;
LineWidth = 0.25;
% f_min_1 = 1; f_max_1 = 150; YTickRange = [1 50 100 150]; XTickRange = -1:3; freqType = 'high';
f_min_1 = 1; f_max_1 = 30; YTickRange = [1 10 20 30]; XTickRange = -1:3; XTickLabel = [-1:3]; freqType = 'low';
method = 'sem'; % 'bootci', 'sem'
valence_label_1 = {'1', '2', '3', '4', '5'};
LineWidth1 = 0.5; LineWidth2 = 0.25;
n_valence = length(valence_label_1);
color_valence = [[235, 80, 49]./255; [245, 204, 196]./255; [200, 200, 200]./255; [62, 164, 201]./255; [34, 21, 117]./255];
freq_ranges = {[13 30]; [60 100]};
freq_labels = {'Beta', 'Gamma'};
for i_range = 1:length(freq_ranges)
    freq_range = freq_ranges{i_range};
    freq_label = freq_labels{i_range};
    for i_reg = 1:n_reg
        close all;
        f = figure('units', 'centimeters', 'position', [5 5 3 2.5], 'Color', [1, 1, 1]); hold on;
        h = cell(n_valence, 1);
        pow_avg_valence = cell(n_valence, 1);
        for i_val = 1:n_valence
            pow = squeeze(pow_valence{i_reg, i_val});
            time_range = [-1 3];
            idxTime = find(freq_tmp.time >= time_range(1) & freq_tmp.time <= time_range(2));
            idxFreq = find(freq_tmp.freq >= freq_range(1) & freq_tmp.freq <= freq_range(2));
            pow_avg = squeeze(mean(pow(:, idxFreq, idxTime), 2)); % trialRepeat_time
            h{i_val} = plot(freq_tmp.time(idxTime), nanmean(pow_avg, 1), 'Color', color_valence(i_val, :), 'LineWidth', 0.5); hold on;
            if strcmp(method, 'bootci')
                ci = bootci(1000, @(x) nanmean(x), pow_avg);
                fill([freq_tmp.time(idxTime) fliplr(freq_tmp.time(idxTime))], [ci(2,:) fliplr(ci(1,:))], color_valence(i_val, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            elseif strcmp(method, 'sem')
                boundedline(freq_tmp.time(idxTime), mean(pow_avg,1,'omitnan'), std(pow_avg,0, 1,'omitnan')./sqrt(size(pow_avg, 2)), ...
                    'cmap', color_valence(i_val, :),'alpha');
            end
            pow_avg_valence{i_val} = pow_avg;
        end
        xlabel('Time (s)')
        ylabel([freq_label ' power']);
        set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth2, 'XColor', 'k', 'YColor', 'k', ...
            'XLim', [-1 3], 'XTick', XTickRange, 'XTickLabel', XTickLabel,...
            'YLim', [-0.4 0.3], 'YTick', [-0.4:0.2:0.2], 'YTickLabel', [-0.4:0.2:0.2], ...
            'tickdir', 'out', 'ticklength', [0.02, 0.02]); box off;


        %--- save data for ART ANOVA in R
        data1 = pow_avg_valence{1}; 
        data2 = pow_avg_valence{2}; 
        data3 = pow_avg_valence{3}; 
        data4 = pow_avg_valence{4}; 
        data5 = pow_avg_valence{5}; 
        data1 = data1(~any(isnan(data1), 2), :);
        data2 = data2(~any(isnan(data2), 2), :);
        data3 = data3(~any(isnan(data3), 2), :);
        data4 = data4(~any(isnan(data4), 2), :);
        data5 = data5(~any(isnan(data5), 2), :);
        time = freq_tmp.time(idxTime); % 1x 21 timepoints
        save([dataName(1:4) '\' dataName '_TFR_' freq_label '_arg_' regions{i_reg} '.mat'], "data1", "data2", "data3", "data4", "data5", "time");

        %--- Trial-level permutation test, time resolved
        offset = [0.01:0.01:0.10]; lindWidth = 1; time = freq_tmp.time(idxTime); % 1x 21 timepoints
        color = repmat([1 1 1], 10, 1);
        lab = yticklabels;
        labV = zeros(size(lab, 1),1);
        for iLab = 1:length(labV)
            labV(iLab) = str2num(lab(iLab, :));
        end
        data_list = {data1, data2, data3, data4, data5};  % cell array of 3 datasets
        [pairs, sig_timepoints_all] = timeResolved_permutation_multiple_conditions(data_list, time, labV, lindWidth, offset, color, 1);

        if max(sig_timepoints_all, [], 'all') > 0
            fprintf('found significance in %s by avging %s band!\n', regions{i_reg}, freq_label)
        end
        
        outputFigure(fullfile([path '\fig3\']), [dataName '_' regions{i_reg} '_' freq_label '_avged']);
    end

end


%% Fig. 3c, beta burst detection
clear -path; close all;
dataName = 'fig3c';
load([path '\fig3\' dataName '.mat' ])
f = figure('units', 'centimeters', 'position', [5, 10, 12, 8], 'Color', [1, 1, 1]);
fontNumbleLable = 6; fontNumbleTick = 6;
LineWidth1 = 0.25; LineWidth2 = 0.5; 
%----------------------------plot data4Ripples
ax1 = subplot(5, 1, 1); % Proprecessed LFP
plot(data4Ripples.time{1}, data4Ripples.trial{1}, '-', 'Color', [0, 0, 0], 'LineWidth', LineWidth1); hold on;
for i_ripple = 1:length(ripples.ripples)
    logIdx = (data4Ripples.time{1} >= ripples.ripples(i_ripple).startTime) & (data4Ripples.time{1} <= ripples.ripples(i_ripple).endTime);
    plot(data4Ripples.time{1}(logIdx), data4Ripples.trial{1}(logIdx), '-', 'Color', rgb('darkgreen'), 'LineWidth', LineWidth2)
end
ylabel({'Macroelec-', 'trode LFP'}, 'FontSize', fontNumbleLable);
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth1, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', [-50 0 50], 'YTickLabel', [-50 0 50], 'XTick', -1:1:3, 'XTickLabel', [], ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
pos1 = get(ax1, 'Position');
subplot(5, 1, 2); % Beta band LFP
plot(data4Ripples.time{1}, data4Ripples.dataB{1}, '-', 'Color', [0, 0, 0], 'LineWidth', LineWidth1); hold on;
for i_ripple = 1:length(ripples.ripples)
    logIdx = (data4Ripples.time{1} >= ripples.ripples(i_ripple).startTime) & (data4Ripples.time{1} <= ripples.ripples(i_ripple).endTime);
    plot(data4Ripples.time{1}(logIdx), data4Ripples.dataB{1}(logIdx), '-', 'Color', rgb('darkgreen'), 'LineWidth', LineWidth2)
end
% xlim([-1 3]);
ylabel({'Beta band', 'LFP'}, 'FontSize', fontNumbleLable);
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth1, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', [-2 0 2], 'YTickLabel', [-2 0 2], 'XTick', -1:1:3, 'XTickLabel', [], ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
% h = gca; h.XAxis.Visible = 'off';
subplot(5, 1, 3); % Beta band envelope
plot(data4Ripples.time{1}, data4Ripples.dataEB{1}, '-', 'Color', [0, 0, 0], 'LineWidth', LineWidth1); hold on;
for i_ripple = 1:length(ripples.ripples)
    logIdx = (data4Ripples.time{1} >= ripples.ripples(i_ripple).startTime) & (data4Ripples.time{1} <= ripples.ripples(i_ripple).endTime);
    plot(data4Ripples.time{1}(logIdx), data4Ripples.dataEB{1}(logIdx), '-', 'Color', rgb('darkgreen'), 'LineWidth', LineWidth2)
end
% xlim([-1 3]);
ylabel({'Beta band', 'envelope'}, 'FontSize', fontNumbleLable);
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth1, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', [0 2 4], 'YTickLabel', [0 2 4], 'XTick', -1:1:3, 'XTickLabel', [], ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
%--- thresholds
center  = mean(data4Ripples.dataEB{1}, 2, 'omitnan');
spread  = std(data4Ripples.dataEB{1}, [], 2, 'omitnan');
threshMin       = center + spread .* param.ripple.threshMinFac; % minimum threshold
threshMax       = center + spread .* param.ripple.threshMaxFac; % maximum threshold
threshPeakMin   = center + spread .* param.ripple.threshPeakMinFac; % minimum peak threshold
yline(threshMin, '--', 'Color', [1 0 0] ,'LineWidth', 0.5);
text(-0.9, threshPeakMin, 'Minimum threshold', 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'FontSize', 6);
% h = gca; h.XAxis.Visible = 'off';
ax4 = subplot(5, 1, 4); % beta band power spectrogram
cfg = [];
cfg.toi = 'all';
cfg.maskstyle = 'saturation';
cfg.parameter = 'powspctrm';
cfg.figure = 'gca';
cfg.title = ' ';
ft_singleplotTFR(cfg, ripplesDetected.tfData);
ylabel({'Power', 'spectrogram'}, 'FontSize', fontNumbleLable);
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth1, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', [15 20 30], 'YTickLabel', [15 20 30], 'XTick', -1:1:3, 'XTickLabel', [], ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
pos4 = get(ax4, 'Position');
set(ax4, 'Position', [pos1(1) pos4(2) pos1(3) pos4(4)]);
colorbarHandle = colorbar; colorbarHandle.Color = 'white'; cb.LineWidth = 0.25; pos = get(colorbarHandle, 'Position');
pos(1) = pos(1) - 0.1; pos(2) = pos(2) + 0.02;
pos(3) = pos(3) * 0.3; pos(4) = pos(4) * 0.7;
set(colorbarHandle, 'Position', pos);
subplot(5, 1, 5); % gamma band LFP
plot(data4RipplesHF.time{1}, data4RipplesHF.dataB{1}, '-', 'Color', [0, 0, 0], 'LineWidth', LineWidth1); hold on;
for i_ripple = 1:length(ripples.ripples)
    logIdx = (data4RipplesHF.time{1} >= ripples.ripples(i_ripple).startTime) & (data4RipplesHF.time{1} <= ripples.ripples(i_ripple).endTime);
    plot(data4RipplesHF.time{1}(logIdx), data4RipplesHF.dataB{1}(logIdx), '-', 'Color', rgb('darkgreen'), 'LineWidth', LineWidth2)
end
% xlim([-1 3]);
xlabel('Time (s)')
ylabel({'Gamma band', 'LFP'}, 'FontSize', fontNumbleLable);
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth1, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', [-5 0 5], 'YTickLabel', [-5 0 5], 'XTick', -1:1:3, 'XTickLabel', {'-1', '0', '1', '2', '3'}, ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;

figName = dataName;
outputFigure(fullfile([path figName(1:4)]), figName);


%% fig. 3d, characteristics in HPC/AMY
dataName_all = {'fig3d_AMY', 'fig3d_HPC'};
for i_fig = 1:length(dataName_all)
    clear -path; close all;
    dataName = dataName_all{i_fig};

    tickRange_all = {[0.2 0.4], [0.05 0.5], [15 30]};
    xLim_all = {[0.1 0.45], [0 0.55], [13 32]};
    load([path '\fig3_update\' dataName '.mat' ])
    valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
    data_all = {duration_group_valence, rate_group_valence, freq_group_valence};
    charLabel_all = {'Burst duration (s)', 'Burst rate (Hz)', 'Burst frequency (Hz)'};
    charLabel_save = {'Duration (s)', 'Rate (Hz)', 'Frequency (Hz)'};
    for i = 1:length(data_all)
        close all;
        data = data_all{i}; charLabel = charLabel_all{i}; tickRange = tickRange_all{i}; xLim = xLim_all{i};
        figName = [dataName '_' charLabel_save{i}];
        f = figure('units', 'centimeters', 'position', [5, 10, 2, 3], 'Color', [1, 1, 1]);
        axes('units', 'centimeters', 'position', [0.5, .9, 1.5, 2]); hold on;
        ax = gca;
        pos = ax.Position;
        fontNumbleLable = 6; fontNumbleTick = 6;
        LineWidth = 0.25;
        visual_type = 'h';
        [D color] = deal(cell(1, n_valence));
        for i_valence = 1:n_valence
            data_tmp = data{i_valence} ;
            D{i_valence} = data_tmp';
            color{i_valence} = color_valence(i_valence, :);
            if i==3
                disp(['burst frequency for valence ', num2str(i_valence) ':']);
                freq_mean = mean(data_tmp)
                freq_std = std(data_tmp)
            end
        end

        makeScatterPlot_switch(D,0,color, 1);
        xlabel(charLabel, 'FontSize', fontNumbleLable);
        ylabel('');
        set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
            'YLim', [0.5 5.5], 'YTick', 1:5, 'YTickLabel', valence_label_1, ...
            'XLim', xLim, 'XTick', tickRange, 'XTickLabel', tickRange, ...
            'tickdir', 'out', 'ticklength', [0.02, 0.02]); box off;

        %--- save data for ART ANOVA in R
        all_data = [];
        all_conditions = [];
        for i = 1:length(D)
            n = numel(D{i});
            all_data = [all_data; D{i}'];
            all_conditions = [all_conditions; repmat(i, n, 1)];
        end
        % Save to .mat file
        save([figName(1:4) '\' figName '.mat'], 'all_data', 'all_conditions');

        %--- Trial-level permutation test
        data_permute = cell2mat(D');
        num_permutations = 1000;
        pvals_fdr = nonTimeResolved_permutation_multiple_conditions(data_permute, num_permutations);
        sig_timepoints = find(pvals_fdr < 0.05);
        disp([figName ', pairs of significance: ']);
        pairs = nchoosek(1:size(data_permute, 1), 2);
        disp(num2str(pairs(sig_timepoints, :)));
        disp('p values:');
        disp(num2str(pvals_fdr(sig_timepoints)));
        outputFigure(fullfile([path 'fig3_update\onset']), figName);
    end

end


function outputFigure(outputpath, figName)
if ~exist(outputpath,'dir')
    mkdir(outputpath); 
end
pdf_file  = fullfile(outputpath, [figName '.pdf']);
if exist(pdf_file, 'file') == 2
    delete(pdf_file);
end
figureHandles = findall(groot, 'Type', 'figure');
figureNumbers = arrayfun(@(h) h.Number, figureHandles);
[~, ind] = sort(figureNumbers);
for i = 1:length(figureNumbers)
    exportgraphics(figureHandles(ind(i)), pdf_file, 'ContentType', 'vector', 'Append', true);
end
end