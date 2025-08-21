%% Fig. 4a, LFPs during beta bursts
dataName_all = {'fig4a_AMY', 'fig4a_HPC'};
method = 'bootci'; % 'bootci', 'sem'
for i_fig = 1:length(dataName_all)
    clear -path; close all;
    dataName = dataName_all{i_fig};
    load([path '\fig4\' dataName '.mat' ]) 
    fontNumbleLable = 6; fontNumbleTick = 6;
    LineWidth1 = 0.5; LineWidth2 = 0.25;
    XTickRange = [-0.1 0 0.2 0.4]; XTickLabel = [-0.1 0 0.2 0.4];

    dataOnRipplesTmp_group_valence = dataOnRipplesHF_group_valence;
    dataOnRipplesTmp_group_valence_surro = dataOnRipplesHF_group_valence_surro;
    dataOnRipplesTmp_group = dataOnRipplesHF_group;
    dataOnRipplesTmp_group_surro = dataOnRipplesHF_group_surro;
    freqTypeLable = 'HG';

    %--- gamma amplitude averaged across all valence conditions
    f = figure('units', 'centimeters', 'position', [5 5 5 7], 'Color', [1, 1, 1]);
    ax1 = axes('Units', 'centimeters', 'Position', [1.2, 5.1, 2.912, 1.5]);
    h2 = plot(burstWindow, mean(dataOnRipplesTmp_group, 1), 'Color', [127 0 255]/255, 'LineWidth', LineWidth1 ); box off; hold on;
    h3 = plot(burstWindow, mean(dataOnRipplesTmp_group_surro, 1), 'Color', [25 0 51]/255, 'LineWidth', LineWidth1 ); hold on;
    if strcmp(method, 'bootci')
        ci = bootci(1000, @(x) mean(x), dataOnRipplesTmp_group);
        fill([burstWindow fliplr(burstWindow)], [ci(2,:) fliplr(ci(1,:))], [127 0 255]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        ci = bootci(1000, @(x) mean(x), dataOnRipplesTmp_group_surro);
        fill([burstWindow fliplr(burstWindow)], [ci(2,:) fliplr(ci(1,:))], [25 0 51]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    elseif strcmp(method, 'sem')
        boundedline(burstWindow, mean(dataOnRipplesTmp_group,1,'omitnan'), std(dataOnRipplesTmp_group,0, 1,'omitnan')./sqrt(size(dataOnRipplesTmp_group, 2)), ...
            'cmap', [127 0 255]/255,'alpha');
        boundedline(burstWindow, mean(dataOnRipplesTmp_group_surro,1,'omitnan'), std(dataOnRipplesTmp_group_surro,0, 1,'omitnan')./sqrt(size(dataOnRipplesTmp_group_surro, 2)), ...
            'cmap', [25 0 51]/255,'alpha');
    end
    ylabel([freqTypeLable ' amplitude']);
    pos1 = get(ax1, 'Position');
    %-- permutation test
    lab = yticklabels;
    labV = zeros(length(yticklabels),1);
    for iLab = 1:length(yticklabels)
        labV(iLab) = str2num(lab{iLab});
    end
    permutation_test_equalSamples(dataOnRipplesTmp_group, dataOnRipplesTmp_group_surro, burstWindow, labV, 1.0, 0.05);
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth2, 'XColor', 'k', 'YColor', 'k', ...
        'YTick', [0 0.5 1], 'YTickLabel', [0 0.5 1],...
        'XLim', [XTickRange(1) XTickRange(end)], 'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
    ax2 = axes('Units', 'centimeters', 'Position', [1.2, 2.9, 3.6, 1.5]);
    imagesc(burstWindow, 1:size(dataOnRipplesTmp_group, 1), dataOnRipplesTmp_group);
    axis xy;
    colorbar;
    climRange = [-1 1];
    set(gca,'YDir','normal', 'CLim', climRange)
    % xlabel('Time (s)');
    ylabel({'Channel number', 'real bursts'})
    pos2 = get(ax2, 'Position');
    set(ax2, 'Position', [pos1(1) pos2(2) pos1(3) pos2(4)]);
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth2, 'XColor', 'k', 'YColor', 'k', ...
        'YTick', [1 15 30 45], 'YTickLabel', [1 15 30 45],'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
    ax3 = axes('Units', 'centimeters', 'Position', [1.2, 0.8, 2.5, 1.5]);
    imagesc(burstWindow, 1:size(dataOnRipplesTmp_group_surro, 1), dataOnRipplesTmp_group_surro);
    axis xy;
    colorbar;
    set(gca,'YDir','normal', 'CLim', climRange)
    xlabel('Time (s)');
    ylabel({'Channel number', 'surrogate bursts'})
    pos3 = get(ax3, 'Position');
    set(ax3, 'Position', [pos1(1) pos3(2) pos1(3) pos3(4)]);
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth2, 'XColor', 'k', 'YColor', 'k', ...
        'YTick', [1 15 30 45], 'YTickLabel', [1 15 30 45],'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;

    figName = dataName;
    outputFigure(fullfile([path figName(1:4)]), [figName '_realVsSurro']);

    %---beta burst-triggerd gamma amplitude on valence condtitions with stats
    dataOnRipplesTmp_group_1 = dataOnRipplesTmp_group_valence{1};
    dataOnRipplesTmp_group_2=[];
    for i_valence = 2:4
        dataOnRipplesTmp_group_2 = cat(1, dataOnRipplesTmp_group_2, dataOnRipplesTmp_group_valence{i_valence});
    end
    dataOnRipplesTmp_group_3 = dataOnRipplesTmp_group_valence{5};
    close all; f = figure('units', 'centimeters', 'position', [5 5 4 3], 'Color', [1, 1, 1]); hold on;
    % axes('Units', 'centimeters', 'Position', [0.5, 0.5, 2.5, 1.5]);
    h1 = plot(burstWindow, mean(dataOnRipplesTmp_group_1,1,'omitnan'), 'Color', color_valence(1, :), 'LineWidth', LineWidth1 ); hold on;
    h2 = plot(burstWindow, mean(dataOnRipplesTmp_group_2,1,'omitnan'), 'Color', [102 204 0]./255, 'LineWidth', LineWidth1 ); hold on;
    h3 = plot(burstWindow, mean(dataOnRipplesTmp_group_3,1,'omitnan'), 'Color', color_valence(5, :), 'LineWidth', LineWidth1 ); hold on;
    if strcmp(method, 'bootci')
        ci = bootci(1000, @(x) mean(x,'omitnan'), dataOnRipplesTmp_group_1);
        fill([burstWindow fliplr(burstWindow)], [ci(2,:) fliplr(ci(1,:))], color_valence(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        ci = bootci(1000, @(x) mean(x,'omitnan'), dataOnRipplesTmp_group_2);
        fill([burstWindow fliplr(burstWindow)], [ci(2,:) fliplr(ci(1,:))], [102 204 0]./255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        ci = bootci(1000, @(x) mean(x,'omitnan'), dataOnRipplesTmp_group_3);
        fill([burstWindow fliplr(burstWindow)], [ci(2,:) fliplr(ci(1,:))], color_valence(5, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    elseif strcmp(method, 'sem')
        boundedline(burstWindow, mean(dataOnRipplesTmp_group_1,1,'omitnan'), std(dataOnRipplesTmp_group_1,0, 1,'omitnan')./sqrt(size(dataOnRipplesTmp_group_1, 2)), ...
            'cmap', color_valence(1, :),'alpha');
        boundedline(burstWindow, mean(dataOnRipplesTmp_group_2,1,'omitnan'), std(dataOnRipplesTmp_group_2,0, 1,'omitnan')./sqrt(size(dataOnRipplesTmp_group_2, 2)), ...
            'cmap', [102 204 0]./255,'alpha');
        boundedline(burstWindow, mean(dataOnRipplesTmp_group_3,1,'omitnan'), std(dataOnRipplesTmp_group_3,0, 1,'omitnan')./sqrt(size(dataOnRipplesTmp_group_3, 2)), ...
            'cmap', color_valence(5, :),'alpha');
    end
    ylabel([freqTypeLable ' amplitude']);
    xlabel('Time (s)');

    %--- save data for ART ANOVA in R
    data1 = dataOnRipplesTmp_group_1; % 30 trials x 21 timepoints
    data2 = dataOnRipplesTmp_group_2; % 120 trials x 21 timepoints
    data3 = dataOnRipplesTmp_group_3; % 30 trials x 21 timepoints
    time = burstWindow; % 1x 21 timepoints
    save([dataName(1:4) '\' dataName '_condition.mat'], "data1", "data2", "data3", "time");

    %--- Trial-level permutation test, time resolved
    offset = [0.05 0.06 0.07]; lindWidth = 1; color_all = [color_valence(1, :); [0, 0, 0]; color_valence(5, :)];
    lab = yticklabels;
    labV = zeros(length(yticklabels),1);
    for iLab = 1:length(yticklabels)
        labV(iLab) = str2num(lab{iLab});
    end
    data_list = {data1, data2, data3};  % cell array of 3 datasets
    [pairs, sig_timepoints_all] = timeResolved_permutation_multiple_conditions(data_list, time, labV, lindWidth, offset, color_all, 1);


    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth2, 'XColor', 'k', 'YColor', 'k', ...
        'YLim', [0 1], 'YTick', [0 0.5 1], 'YTickLabel', [0 0.5 1],...
        'XLim', [XTickRange(1) XTickRange(end)], 'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;

    figName = dataName;
    outputFigure(fullfile([path figName(1:4)]), [figName '_condition']);

    %--- beta burst-triggerd beta ITPC on valence condtitions with stats
    close all;
    ITPCOnRipplesTmp_group_valence = ITPCOnRipplesLF_group_valence;
    ITPCOnRipplesTmp_group_valence_surro = ITPCOnRipplesLF_group_valence_surro;
    ITPCOnRipplesTmp_group_1 = ITPCOnRipplesTmp_group_valence{1};
    color_1 = color_valence(1, :);
    ITPCOnRipplesTmp_group_2 = [];
    for i_valence = 1:4
        ITPCOnRipplesTmp_group_2 = cat(1, ITPCOnRipplesTmp_group_2, ITPCOnRipplesTmp_group_valence{i_valence});
    end
    color_2 = [102 204 0]./255;
    ITPCOnRipplesTmp_group_3 = ITPCOnRipplesTmp_group_valence{5};
    color_3 = color_valence(5, :);
    legendTxt = {'1', '234', '5'};
    data_group_1 = ITPCOnRipplesTmp_group_1;
    data_group_2 = ITPCOnRipplesTmp_group_2;
    data_group_3 = ITPCOnRipplesTmp_group_3;
    limRangeVal = [-0.3 0.3];
    text_1 = ['ITPC during ' tickLabel ' burst'];
    text_2 = 'ITPC';
    % burstWindow_2 = -0.2:0.02:0.2;
    burstWindow_2 = -0.1:0.02:0.4;
    close all; f = figure('units', 'centimeters', 'position', [5 5 4 3], 'Color', [1, 1, 1]); hold on;
    h1 = plot(burstWindow_2, squeeze(nanmean(data_group_1, [1 2])), 'Color', color_1, 'LineWidth', LineWidth1); hold on;
    h2 = plot(burstWindow_2, squeeze(nanmean(data_group_2, [1 2])), 'Color', color_2, 'LineWidth', LineWidth1); hold on;
    h3 = plot(burstWindow_2, squeeze(nanmean(data_group_3, [1 2])), 'Color', color_3, 'LineWidth', LineWidth1); hold on;
    if strcmp(method, 'bootci')
        ci = bootci(1000, @(x) nanmean(x), squeeze(nanmean(data_group_1, 2)));
        fill([burstWindow_2 fliplr(burstWindow_2)], [ci(2,:) fliplr(ci(1,:))], color_1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        ci = bootci(1000, @(x) nanmean(x), squeeze(nanmean(data_group_2, 2)));
        fill([burstWindow_2 fliplr(burstWindow_2)], [ci(2,:) fliplr(ci(1,:))], color_2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        ci = bootci(1000, @(x) nanmean(x), squeeze(nanmean(data_group_3, 2)));
        fill([burstWindow_2 fliplr(burstWindow_2)], [ci(2,:) fliplr(ci(1,:))], color_3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    elseif strcmp(method, 'sem')
        boundedline(burstWindow_2, mean(squeeze(nanmean(data_group_1, 2)),1,'omitnan'), ...
            std(squeeze(nanmean(data_group_1, 2)),0, 1,'omitnan')./sqrt(size(squeeze(nanmean(data_group_1, 2)), 2)), 'cmap', color_1,'alpha');
        boundedline(burstWindow_2, mean(squeeze(nanmean(data_group_2, 2)),1,'omitnan'), ...
            std(squeeze(nanmean(data_group_2, 2)),0, 1,'omitnan')./sqrt(size(squeeze(nanmean(data_group_2, 2)), 2)), 'cmap', color_2,'alpha');
        boundedline(burstWindow_2, mean(squeeze(nanmean(data_group_3, 2)),1,'omitnan'), ...
            std(squeeze(nanmean(data_group_3, 2)),0, 1,'omitnan')./sqrt(size(squeeze(nanmean(data_group_3, 2)), 2)), 'cmap', color_3,'alpha');
    end
    xlabel('Time (s)');
    ylabel(text_2)
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth2, 'XColor', 'k', 'YColor', 'k', ...
        'XLim', [XTickRange(1) XTickRange(end)], 'XTick', XTickRange, 'XTickLabel', XTickLabel,...
        'YLim', [0.15 0.45], 'YTick', [0.15 0.3 0.45], 'YTickLabel', [0.15 0.3 0.45], ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02]); box off;
    %--- save data for ART ANOVA in R
    data1 = squeeze(nanmean(data_group_1, 2)); % 30 trials x 21 timepoints
    data2 = squeeze(nanmean(data_group_2, 2)); % 120 trials x 21 timepoints
    data3 = squeeze(nanmean(data_group_3, 2)); % 30 trials x 21 timepoints
    time = burstWindow_2; % 1x 21 timepoints
    save([dataName(1:4) '\' dataName '_ITPC.mat'], "data1", "data2", "data3", "time");

    %--- Trial-level permutation test, time resolved
    offset = [0.01 0.02 0.03]; lindWidth = 1; color_all = [color_valence(1, :); [0, 0, 0]; color_valence(5, :)];
    lab = yticklabels;
    labV = zeros(size(lab, 1),1);
    for iLab = 1:length(labV)
        labV(iLab) = str2num(lab(iLab, :));
    end
     data_list = {data1, data2, data3};  % cell array of 3 datasets
    [pairs, sig_timepoints_all] = timeResolved_permutation_multiple_conditions(data_list, time, labV, lindWidth, offset, color_all, 1);

    figName = dataName;
    outputFigure(fullfile([path figName(1:4) ]), [figName '_ITPC']);

end

%% Fig. 4b, FRs during beta bursts
figName_all = {'fig4b1', 'fig4b2', 'fig4c1', 'fig4c2', 'fig4d1', 'fig4d2', 'fig4e1', 'fig4e2'};
ylim = [-0.04 0.15; -0.1 0.22; -0.0 0.12; -0.02 0.22; 0 0.12; -0.1 0.3; 0 0.15; -0.15 0.26];
for i_fig = 1:length(figName_all)
    figName = figName_all{i_fig};
    clear -path; close all;
    fontNumbleLable = 6; fontNumbleTick = 6;
    LineWidth = 0.25;
    YTickRange = [0 0.5 1]; YTickLabel = [0 0.5 1]; XTickRange = [-0.1 0 0.2 0.4]; XTickLabel = [-0.1 0 0.2 0.4];
    fig = openfig([path 'fig4\' figName '.fig']); box off;
    ax = findall(fig, 'type', 'axes');
    for i = 1:length(ax)
        lines = findall(ax(i), 'type', 'line');
        for j = 1:length(lines)
            if contains(lines(j).DisplayName, 'data') || isempty(lines(j).DisplayName)
                lines(j).LineWidth = 1.0;
                if strcmp(figName, 'fig4b1')
                    y = lines(i).YData;
                    lines(j).YData = y - 0.025;
                elseif strcmp(figName, 'fig4e1')
                    y = lines(i).YData;
                    lines(j).YData = y - 0.01;
                end
            else
                lines(j).LineWidth = 0.5;
            end
        end
        ax(i).Position = [0.18 0.21 0.7 0.8];
    end
    set(fig, 'Units', 'centimeters', 'Position', [4 4 4 3], 'Color', [1 1 1]);
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
        'YLim', ylim(i_fig, :), 'YTick', [0 0.1 0.2], 'YTickLabel', [0 0.1 0.2],...
        'XLim', [-0.1 0.4], 'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]);
    title(gca, '');
    ylabel('Z-scored FR')
    xlabel('Time (s)')
    lgd = legend('Location','northeast'); legend('boxoff');
    lgd.FontSize = 6;
    lgd.ItemTokenSize = [6, 18];
    h = findobj(gcf, 'Type', 'Legend');
    set(h, 'Visible', 'off');

    outputFigure(fullfile([path figName(1:4)]), figName(1:6));
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