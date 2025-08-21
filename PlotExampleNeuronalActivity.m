%% Raster part
spike_time_1 = [];
xspikes_ple_1 = []; xspikes_unple_1 = []; xspikes_neu_1 = [];
yspikes_ple_1 = []; yspikes_unple_1 = []; yspikes_neu_1 = [];
[N_all, FR_all] = deal([]);
[xspikes_valence, yspikes_valence, N_valence, FR_valence, FR_valence_smoothed, FR_valence_z] = deal(cell(n_valence, 1));
for i_trial = 1:n_trials

    ind = rawSpikeTime >= trigTime(i_trial)+ trial_duration(1) & rawSpikeTime <= trigTime(i_trial)+trial_duration(2);
    xspikes_1 = repmat(rawSpikeTime(ind)' - trigTime(i_trial), 2, 1);
    yspikes_1 = nan(size(xspikes_1));


    %----------- count spikes for different conditions
    %--- regardless of condition
    N_all = cat(1, N_all, histcounts(xspikes_1(1,:),trialTime,'Normalization','count'));

    %--- for valence: 'Very Pleasant', 'Pleasant', 'Neutral', 'Unpleasant', 'Very Unpleasant'
    for i_valence = 1:length(trials_valence)
        if ismember(i_trial, trials_valence{i_valence})
            j_trial = find(trials_valence{i_valence} == i_trial);
            xspikes_valence{i_valence} = [xspikes_valence{i_valence} xspikes_1];
            yspikes_1(1, :) = j_trial - 1;
            yspikes_1(2, :) = j_trial;
            yspikes_valence{i_valence} = [yspikes_valence{i_valence} yspikes_1];
            N_valence{i_valence} = cat(1, N_valence{i_valence}, histcounts(xspikes_1(1,:),trialTime,'Normalization','count'));
        end
    end

end

%----------- raster plot for different conditions
xOffset = 0.26; yOffset1 = 0.39; yOffset2 = 0.13;
bh1 = 0.53; ih = 0.015; bw = 0.580; bh2 = 0.25;
f = figure('units', 'centimeters', 'position', [5 5 4.0 4.5], 'Color', [1, 1, 1]);
% bh1 = 0.53; ih = 0.015; bw = 0.650; bh2 = 0.25;
% f = figure('units', 'centimeters', 'position', [5 5 3.5 4.5], 'Color', [1, 1, 1]);
% ax1 = axes('Units', 'centimeters', 'Position', [1.2, 5.1, 2.912, 1.5]);
for i_valence = 1:n_valence
    if i_valence > 1
        prob_tmp1 = sum(cellfun(@length, trials_valence(1:i_valence - 1)))/n_trials;
    else
        prob_tmp1 = 0;
    end
    prob_tmp2 = length(trials_valence{i_valence})/n_trials;
    H = subplot('position', [xOffset, yOffset1 + bh1 * prob_tmp1 + ih * (i_valence - 1), bw, bh1 * prob_tmp2]);
    plot(xspikes_valence{i_valence}, yspikes_valence{i_valence}, 'Color', color_valence(i_valence, :), 'LineWidth', 0.25);
    h = gca;
    h.XAxis.Visible = 'off';
    box off;
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
        'XLim', [-1 3], 'XTick', XTickRange, 'XTickLabel', XTickLabel,...
        'YLim', [0, length(trials_valence{i_valence})], 'YTick', [1, length(trials_valence{i_valence})], 'YTickLabel', [1, length(trials_valence{i_valence})], ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
    xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);
    text(H.XLim(2) + 0.2, mean(H.YLim), valence_label_1{i_valence}, 'Color', color_valence(i_valence, :), 'FontSize', 6)
end
h = ylabel({'Trials'}, 'FontSize', fontNumbleLable);
pos = get(h, 'Position');
set(h, 'Position', [pos(1), pos(2)-55, pos(3)]);

%% FR part
FR_all = N_all/step_size;
FR_all_smoothed = smoothdata(FR_all,2,'gaussian', bin_size/step_size);

%----------- baseline correction
% tw_norm = 1:1/step_size;
tw_norm = find(trialTime >= -1 & trialTime <= 0);
baselineType = 'absolute'; % 'absolute', 'relative', 'relchange', 'normchange', 'db', 'zscore', 'no'

FR_all_z = firingRate_normalization(FR_all_smoothed, baselineType, tw_norm);

%------------ FR plot
H0 = subplot('position', [xOffset, yOffset2, bw, bh2]); hold on;
box off
xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);
xlabel('Time (s)', 'FontSize', fontNumbleLable);
ylabel({'Relative FR (Hz)'}, 'FontSize', fontNumbleLable)
for i_valence = 1:n_valence
    FR_valence{i_valence} = N_valence{i_valence}/step_size;
    FR_valence_smoothed{i_valence} = smoothdata(FR_valence{i_valence},2,'gaussian', bin_size/step_size);
    FR_valence_z{i_valence} = firingRate_normalization(FR_valence_smoothed{i_valence}, baselineType, tw_norm);
    boundedline(trialTime(1:end-1), mean(FR_valence_z{i_valence}, 1), std(FR_valence_z{i_valence}, 1)./sqrt(size(FR_valence_z{i_valence}, 1)), 'cmap', color_valence(i_valence, :),'alpha'); hold on
    % plot(trialTime(1:end-1), mean(FR_valence_z{i_valence}, 1), 'Color', color_valence(i_valence, :), 'LineWidth', 0.2);
end
ax = findall(H0, 'type', 'axes');
for i = 1:length(ax)
    lines = findall(ax(i), 'type', 'line');
    for j = 1:length(lines)
        lines(j).LineWidth = 0.25; 
    end
end
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
    'XLim', [-1 3], 'XTick', XTickRange, 'XTickLabel', XTickLabel,...
    'YLim', YLim, 'YTick', YTickRange, 'YTickLabel', YTickLabel, ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);

lab = yticklabels;
labV = zeros(size(lab, 1),1);
for iLab = 1:length(labV)
    labV(iLab) = str2num(lab(iLab, :));
end
minV_val = min(labV);

%% permutation testing
%------ choose the task window and the baseline window
[~, idx_1] = min(abs(trialTime - rp_window(1))); 
[~, idx_2] = min(abs(trialTime - rp_window(2))); 
tw_task = idx_1:idx_2;
[~, idx_1] = min(abs(trialTime - bl_window(1))); 
[~, idx_2] = min(abs(trialTime - bl_window(2))); 
tw_bs = idx_1:idx_2;

%--- within valence: 'Very Pleasant', 'Pleasant', 'Neutral', 'Unpleasant', 'Very Unpleasant'
task_mark = zeros(n_valence, 1);
for i_valence = 1:n_valence

    tmp        = struct();
    tmp.time   = trialTime(1:end-1);
    tmp.label  = {'foo'};
    tmp.dimord = 'rpt_chan_time';
    d1 = tmp;
    d2 = tmp;
    d1.time = trialTime(tw_task);
    d2.time = trialTime(tw_task);
    d1.trial(:,1,:) = FR_valence_z{i_valence}(:, tw_task);
    d2.trial(:,1,:) = repmat(mean(FR_valence_z{i_valence}(:, tw_bs), 2), [1, length(tw_task)]);
    pos_sig   = minV_val + 0.05 + 0.25*(i_valence - 1);
    stat = permutation_test_individual(d1, d2);
    sigline = nan(1,numel(tmp.time));
    idx = nearest(tmp.time,stat.time(1)):nearest(tmp.time,stat.time(end));
    % sigline(idx(stat.mask==1)) = pos_sig;
    sigline(idx(stat.prob<0.05)) = pos_sig;
    hold on
    % yyaxis left
    % figure(f0);
    plot(tmp.time,sigline,'-','color',color_valence(i_valence, :),'LineWidth',1)

end



function FR = firingRate_normalization(FR_pre, baselineType, tw_norm)
switch baselineType
    case 'zscore'
        FR_bl_avg_tmp = repmat(mean(FR_pre(:, tw_norm), 2), [1 size(FR_pre, 2)]);
        FR_bl_std_tmp = repmat(std(FR_pre(:, tw_norm), 1, 2), [1 size(FR_pre, 2)]);
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


function stat = permutation_test_individual(fr1, fr2)
n_trial = size(fr1.trial, 1);
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
stat = ft_timelockstatistics(cfg, fr1, fr2);
end


