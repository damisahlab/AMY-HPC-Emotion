n_valence = 5;
responsive_marks = {'enhanced', 'suppressed'};
for respNeuronIdx = 1:n_valence

    for iRepMode = 1:length(responsive_mode)
        close all;
        f = figure('units', 'centimeters', 'position', [5 5 4 3], 'Color', [1, 1, 1]); hold on;
        box off
        xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);
        xlabel('Time (s)', 'FontSize', 6);
        ylabel('Normalized FR (Hz)', 'FontSize', 6)
        set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
            'XLim', XLim, 'XTick', XTickRange, 'XTickLabel', XTickLabel,...
            'YLim', YLim(iRepMode, :), 'YTick', YTickRange(iRepMode, :), 'YTickLabel', YTickLabel(iRepMode, :), ...
            'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
        % title(['Pref ' num2str(valence_number(respNeuronIdx))], 'FontSize', fontNumbleLable);

        switch responsive_mode{iRepMode}
            case 'excitatory'
                ind_tmp = indices_valence.excitatory{respNeuronIdx};
            case 'inhibitory'
                ind_tmp = indices_valence.inhibitory{respNeuronIdx};
        end

        if ~isempty(ind_tmp)
            FR_valence_all_tmp_1 = cell2mat(FR_valence_all{respNeuronIdx}(ind_tmp)');
            FR_valence_all_tmp_2 = [];
            comparNeuronIdx = setdiff(1:5, respNeuronIdx);
            color_comp = [102 204 0]./255;

            for i_valence = 1:length(comparNeuronIdx)
                FR_valence_all_tmp_2 = cat(1, FR_valence_all_tmp_2, cell2mat(FR_valence_all{comparNeuronIdx(i_valence)}(ind_tmp)'));
            end
            color1 = color_valence(respNeuronIdx, :); color2 = color_comp;
            % legend1 = num2str(valence_number(respNeuronIdx)); legend2 = num2str(valence_number(comparNeuronIdx));
            legend1 = valence_label_1{respNeuronIdx}; legend2 = 'Other';
            plot_FR_groups(FR_valence_all_tmp_1, FR_valence_all_tmp_2, trialTime, color1, color2, legend1, legend2)

            outputFigure(fullfile([path dataName(1:4)]), [dataName '_Pref' num2str(valence_number(respNeuronIdx)) '_' responsive_marks{iRepMode} '_1']);

            if respNeuronIdx == 1 || respNeuronIdx == 5
                close all;
                f = figure('units', 'centimeters', 'position', [5 5 4 3], 'Color', [1, 1, 1]); hold on;
                box off
                xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);
                xlabel('Time (s)', 'FontSize', 6);
                ylabel('Normalized FR (Hz)', 'FontSize', 6)
                set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
                    'XLim', XLim, 'XTick', XTickRange, 'XTickLabel', XTickLabel,...
                    'YLim', YLim(iRepMode, :), 'YTick', YTickRange(iRepMode, :), 'YTickLabel', YTickLabel(iRepMode, :), ...
                    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
                % title(['Pref ' num2str(valence_number(respNeuronIdx))], 'FontSize', fontNumbleLable);
                if respNeuronIdx == 1
                    comparNeuronIdx = 5;
                else
                    comparNeuronIdx = 1;
                end
                color_comp = color_valence(comparNeuronIdx, :);
                FR_valence_all_tmp_2 = [];
                for i_valence = 1:length(comparNeuronIdx)
                    FR_valence_all_tmp_2 = cat(1, FR_valence_all_tmp_2, cell2mat(FR_valence_all{comparNeuronIdx(i_valence)}(ind_tmp)'));
                end
                color1 = color_valence(respNeuronIdx, :); color2 = color_comp;
                % legend1 = num2str(valence_number(respNeuronIdx)); legend2 = num2str(valence_number(comparNeuronIdx));
                legend1 = valence_label_1{respNeuronIdx}; legend2 = valence_label_1{comparNeuronIdx};
                plot_FR_groups(FR_valence_all_tmp_1, FR_valence_all_tmp_2, trialTime, color1, color2, legend1, legend2);

                outputFigure(fullfile([path dataName(1:4)]), [dataName '_Pref' num2str(valence_number(respNeuronIdx)) '_' responsive_marks{iRepMode} '_2']);
            end

        end 

    end

end


%% Supportive functions
function plot_FR_groups(FR_valence_all_tmp_1, FR_valence_all_tmp_2, trialTime, color1, color2, legend1, legend2)
tw_norm = find(trialTime >= -0.5 & trialTime <= 0);
baselineType = 'absolute'; % 'absolute', 'relative', 'relchange', 'normchange', 'db', 'zscore', 'no'
FR_valence_all_tmp_1 = firingRate_normalization(FR_valence_all_tmp_1, baselineType, tw_norm);
FR_valence_all_tmp_2 = firingRate_normalization(FR_valence_all_tmp_2, baselineType, tw_norm);

ci = bootci(1000, @(x) mean(x,'omitnan'), FR_valence_all_tmp_1);
fill([trialTime(1:end-1) fliplr(trialTime(1:end-1))], [ci(2,:) fliplr(ci(1,:))], color1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p1 = plot(trialTime(1:end-1), nanmean(FR_valence_all_tmp_1, 1), 'Color', color1, 'LineWidth',0.5);
ci = bootci(1000, @(x) mean(x,'omitnan'), FR_valence_all_tmp_2);
fill([trialTime(1:end-1) fliplr(trialTime(1:end-1))], [ci(2,:) fliplr(ci(1,:))], color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(trialTime(1:end-1), nanmean(FR_valence_all_tmp_2, 1), 'Color', color2, 'LineWidth',0.5);

%-- permutation test
lab = yticklabels;
labV = zeros(size(lab, 1),1);
for iLab = 1:length(labV)
    labV(iLab) = str2num(lab(iLab, :));
end
minV_val = min(labV);
permutation_test_unequalSamples(FR_valence_all_tmp_1, FR_valence_all_tmp_2, trialTime(1:end-1), labV);

lgd = legend([p1, p2], {legend1, legend2});
legend('boxoff');
lgd.FontSize = 6;
% lgd.Position(1) = lgd.Position(1) + 0.25;
lgd.ItemTokenSize = [6, 18];

end


function permutation_test_unequalSamples(data1, data2, time, labV)
n_iter = 5;
p_iter = zeros(n_iter, length(time));
for i_iter = 1:n_iter
    if size(data1, 1) < size(data2, 1)
        ind_1 = 1:size(data1, 1);
        ind_2 = randperm(size(data2, 1), size(data1, 1));
    else
        ind_2 = 1:size(data2, 1);
        ind_1 = randperm(size(data1, 1), size(data2, 1));
    end
    tmp        = struct();
    tmp.time   = time;
    tmp.label  = {'foo'};
    tmp.dimord = 'rpt_chan_time'; % catDim_1_time
    d1 = tmp;
    d2 = tmp;
    d1.trial(:,1,:) = data1(ind_1, :);
    d2.trial(:,1,:) = data2(ind_2, :);
    % color_sig = 'k';
    % pos_sig = labV(1) + 0.15*(labV(2) - labV(1));
    n_trial = size(d1.trial, 1);
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    % cfg.latency             = trial_duration;
    cfg.latency             = 'all';
    cfg.statistic           = 'depsamplesT';
    cfg.alpha               = 0.05;
    cfg.clusteralpha        = 0.1;
    cfg.correcttail         = 'alpha';
    cfg.tail                = 0;
    cfg.clustertail         = cfg.tail;
    cfg.numrandomization    = 1000;
    cfg.neighbours          = [];
    design = zeros(2, 2*n_trial);
    design(1, :) = [ones(1, n_trial) 2*ones(1, n_trial)];
    design(2, :) = [1:n_trial 1:n_trial];
    cfg.design              = design;
    cfg.ivar = 1;
    cfg.uvar = 2;
    stat = ft_timelockstatistics(cfg, d1, d2);
    p_iter(i_iter, :) = stat.prob;
end
prob_avg = mean(p_iter,1);
sigline = nan(1,numel(tmp.time));
idx = nearest(tmp.time,stat.time(1)):nearest(tmp.time,stat.time(end));
sigline(idx(prob_avg < 0.05)) = labV(1) + 0.05;
hold on
plot(tmp.time,sigline, '-k', 'LineWidth',1)
end


function FR = firingRate_normalization(FR_pre, baselineType, tw_norm)
switch baselineType
    case 'zscore'
        % % FR_bl_avg = mean(FR_pre(:, tw_norm));
        % % FR_bl_std = std(FR_pre(:, tw_norm));
        % % FR = (FR_pre - FR_bl_avg) / FR_bl_std;
        % epsilon = 1e-6;
        FR_bl_avg_tmp = repmat(mean(FR_pre(:, tw_norm), 2), [1 size(FR_pre, 2)]);
        FR_bl_std_tmp = repmat(std(FR_pre(:, tw_norm), 1, 2), [1 size(FR_pre, 2)]);
        % FR_bl_std_tmp = max(FR_bl_std_tmp, epsilon);
        FR = (FR_pre - FR_bl_avg_tmp)./ FR_bl_std_tmp;
    case 'absolute'
        FR_bl_avg_tmp = repmat(mean(FR_pre(:, tw_norm), 2), [1 size(FR_pre, 2)]);
        FR = FR_pre - FR_bl_avg_tmp;
    case 'normchange'
        FR_bl_min = repmat(min(FR_pre(:, tw_norm), [], 2), [1 size(FR_pre, 2)]);
        FR_bl_max = repmat(max(FR_pre(:, tw_norm), [], 2), [1 size(FR_pre, 2)]);
        FR = (FR_pre - FR_bl_min) ./ (FR_bl_max - FR_bl_min);
    case 'no'
        FR = FR_pre;
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

