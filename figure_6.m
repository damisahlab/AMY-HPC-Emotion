%% Fig. 6a, representative single-trial of PAC
clear -path; close all;
dataName = 'fig6a';
load([path '\fig6\' dataName '.mat' ])
f = figure('units', 'centimeters', 'position', [5, 5, 5, 4], 'Color', [1, 1, 1]);
hold on;
box off
i_trial = 1
beta = squeeze(beta_all.data(:, :, i_trial));
gamma = squeeze(gamma_all.data(:, :, i_trial));
% Compute theta phase using Hilbert transform
theta_phase = angle(hilbert(beta));
% Compute gamma amplitude envelope
gamma_amp = abs(hilbert(gamma));
% Compute PAC (Phase-Amplitude Coupling)
pac_signal = gamma_amp .* cos(theta_phase); % Modulated signal
% Plot raw LFP
linewidth1 = 0.5; FontSize = 6;
subplot(4,1,1);
plot(t, lfp(:, i_trial)', 'k', 'LineWidth', linewidth1);
xlim([t(1), t(1000)]); % show only a segment
axis tight;
set(gca, 'XColor', 'none', 'YColor', 'none'); % Keep y-label only
set(gca, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'Box', 'off');
ylabel('Raw', 'Color', 'k', 'FontSize', FontSize);
% Plot beta
subplot(4,1,2);
plot(t, beta, 'k', 'LineWidth', linewidth1);
xlim([t(1), t(1000)]);
axis tight;
set(gca, 'XColor', 'none', 'YColor', 'none'); % Keep y-label only
set(gca, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'Box', 'off');
ylabel('Beta', 'Color', 'k', 'FontSize', FontSize);
% axis off;
% Plot gamma
subplot(4,1,3);
plot(t, gamma, 'k', 'LineWidth', linewidth1);
xlim([t(1), t(1000)]);
axis tight;
set(gca, 'XColor', 'none', 'YColor', 'none'); % Keep y-label only
set(gca, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'Box', 'off');
ylabel('Gamma', 'Color', 'k', 'FontSize', FontSize);
% Plot PAC with colored amplitude
subplot(4,1,4);
colors = hot(256);
amp_scaled = rescale(gamma_amp, 1, 256); % scale amplitude for color map
hold on;
for i = 1:999
    color_idx = round(amp_scaled(i));
    plot(t(i:i+1), beta(i:i+1), 'Color', colors(color_idx,:), 'LineWidth', linewidth1);
end
hold off;
xlim([t(1), t(1000)]);
axis tight;
set(gca, 'XColor', 'none', 'YColor', 'none'); % Keep y-label only
set(gca, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', [], 'Box', 'off');
ylabel('PAC', 'Color', 'k', 'FontSize', FontSize);

figName = dataName;
outputFigure(fullfile([path figName(1:4)]), figName);

%% Fig. 6b, phase
clear -path; close all;
dataName = 'fig5b_HPC';
% dataName = 'fig5b_AMY';
load([path '\fig5_update\' dataName '.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
fontNumbleTick = 6; LineWidth = 0.25;
%-- phase distribution
close all;
for i_reg_a = i_start_a:i_end_a
    for i_reg_p = i_start_p:i_end_p

        f = figure('units', 'centimeters', 'position', [5, 5, 14, 3], 'Color', [1, 1, 1]); axis off; hold on;
        all_theta = [];
        all_group_labels = [];
        for i_valence = 1:n_valence

            %---polar plot, method 1
            phase2power_LF = squeeze(mean(phase2power_LF_valence{i_valence, i_reg_p, i_reg_a}, 1));
            % mean_amplitude_per_phase: average over frequencies
            mean_amplitude_per_phase = mean(phase2power_LF, 1); % 1 x 18
            std_amplitude_per_phase = mean(phase2power_LF, 1); % 1 x 18
            amp_max = max(mean_amplitude_per_phase);
            amp_min = min(mean_amplitude_per_phase);
            mean_amplitude_per_phase1 = (mean_amplitude_per_phase - amp_min) / (amp_max - amp_min);
            % Prepare phase bins
            phase_edges_deg = linspace(-180, 180, n_bins+1); % Bin edges
            phase_centers_deg = (phase_edges_deg(1:end-1) + phase_edges_deg(2:end))/2;
            phase_centers_rad = deg2rad(phase_centers_deg);
            % --- Create correct pseudo_theta ---
            scaling_factor = 1000;  % can adjust for smoothing (recommended 300~1000)
            max_value = max(mean_amplitude_per_phase1);
            pseudo_theta = []; reps_all = nan(n_bins, 1);
            for i = 1:n_bins
                reps = max(round(mean_amplitude_per_phase1(i) / max_value * scaling_factor), 1); % proportionally scale
                pseudo_theta = [pseudo_theta; repmat(phase_centers_rad(i), reps, 1)];
                reps_all(i) = reps;
            end

            all_theta = [all_theta; pseudo_theta];
            all_group_labels = [all_group_labels; repmat(i_valence, size(pseudo_theta))];


            % --- PLOT ---
            pos_width = 1 / n_valence * 0.6;
            pos_left = (i_valence - 1) * pos_width;
            ax_pos = [0.035 + (i_valence - 1) * 1 / n_valence * 1.02, 0.05, pos_width, 0.7];  % [left, bottom, width, height] in normalized units
            pax = polaraxes('Position', ax_pos);  % Place polaraxes manually
            hold(pax, 'on')
            hold on
            h = polarhistogram(pax, pseudo_theta, deg2rad(phase_edges_deg), ...
                'Normalization', 'probability', ...
                'FaceColor', color_valence(i_valence, :), 'FaceAlpha', 0.6, 'EdgeColor', color_valence(i_valence, :), 'LineWidth', 0.5);
            % Beautify
            pax.ThetaZeroLocation = 'top';
            pax.ThetaDir = 'clockwise';
            pax.ThetaTick = [0 90 180 270];
            pax.ThetaTickLabel = {'0', '\pi/2', '\pi', '3\pi/2'};
            pax.FontSize = 6;
            pax.RTick = [0 0.04 0.08];
            pax.ThetaColor = [0 0 0];
            pax.RColor = [0 0 0];
            pax.LineWidth = 0.25;
            title(valence_label_1{i_valence}, 'FontSize', fontNumbleTick);
            [mu, ul, ll] = circ_mean(pseudo_theta);
            [pval, z] = circ_rtest(pseudo_theta);
            fprintf('p value for valence %s : %f \n', num2str(i_valence), pval)
            polarplot([mu mu], [0, mean(h.Values)], 'k', 'LineWidth', 1)
            hold off


        end


        % Run the test
        [p_ww, table_ww] = circ_kuipertest(all_theta, all_group_labels);
        disp(['p-value = ', num2str(p_ww)]);


        %--- pairwise Mardia-Watson-Wheeler tests with FDR correction
        % Unique group labels
        group_ids = unique(all_group_labels);
        num_groups = numel(group_ids);
        % Store results
        pvals = zeros(nchoosek(num_groups, 2), 1);
        pairs = zeros(nchoosek(num_groups, 2), 2);
        pair_idx = 1;
        % Loop over all pairwise combinations
        for i = 1:num_groups
            for j = i+1:num_groups
                % Extract data for each group
                theta1 = all_theta(all_group_labels == group_ids(i));
                theta2 = all_theta(all_group_labels == group_ids(j));
                % Mardia-Watson-Wheeler test
                % p = circ_mtest(theta1, theta2);  % You can use circ_kuipertest or circ_cmtest too
                p = circ_kuipertest(theta1, theta2);
                % res = circ_cmtest({theta1, theta2});
                pvals(pair_idx) = p;
                pairs(pair_idx, :) = [group_ids(i), group_ids(j)];
                pair_idx = pair_idx + 1;
            end
        end
        % FDR correction (Benjamini-Hochberg)
        [~, ~, ~, adj_pvals] = fdr_bh(pvals);
        % Display results
        fprintf('Pairwise Mardia-Watson-Wheeler Tests with FDR Correction:\n');
        for k = 1:length(pvals)
            fprintf('Group %d vs Group %d: raw p = %.4f, FDR p = %.4f\n', ...
                pairs(k,1), pairs(k,2), pvals(k), adj_pvals(k));
        end

    end
end

outputFigure(fullfile([path dataName(1:4)]), [dataName '_phase']);

%% Fig. 6c, PAC comodulograms 
clear -path; close all;
dataName = 'fig6b_HPC';
% dataName = 'fig6b_AMY';
fontNumbleTick = 6; LineWidth = 0.25;

YTickRange = [60 100 140]; XTickRange = [15 20 25 30]; XTickLabel = {'', '20', '', '30'};
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
for i_reg_a = i_start_a:i_end_a
    for i_reg_p = i_start_p:i_end_p

        f = figure('units', 'centimeters', 'position', [5, 5, 20, 2], 'Color', [1, 1, 1]); hold on;
        for i_valence = 1:n_valence
            phase2power_HFLF = mean(phase2power_valence{i_valence, i_reg_p, i_reg_a}, 1)';
            comdlgrm_tmp = squeeze(mean(comdlgrm_valence{i_valence, i_reg_p, i_reg_a}, 1));
            phase2power_LF = squeeze(mean(phase2power_LF_valence{i_valence, i_reg_p, i_reg_a}, 1));

            subplot(1,n_valence,i_valence)

            contourf(LF_steps, HF_steps,comdlgrm_tmp',500,'linecolor','none')
            xlabel('Freq. for phase (Hz)')
            if i_valence == 1
            ylabel({'Freq. for ampli-' 'tude (Hz)'})
            end
            if i_valence < 5
                colorbar off
            else
                c = colorbar; 
            end
            set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
            'YTick', YTickRange, 'YTickLabel', YTickRange,'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
            'tickdir', 'out', 'ticklength', [0.01, 0.01], 'clim',[0 1.5e-4]); 
            % box off;
            
            axis square
            title(valence_label_1{i_valence}, FontSize=fontNumbleTick);
        end

    end

end

outputFigure(fullfile([path dataName(1:4)]), dataName);


%% Fig. 6d, stat
clear -path; close all;
dataName = 'fig5b_HPC';
% dataName = 'fig5b_AMY';
load([path '\fig5_update\' dataName '.mat' ])
fontNumbleTick = 6; LineWidth = 0.25;
YTickRange = [0.00002, 0.0001, 0.0002]; XTickRange = [15 20 25 30]; XTickLabel = {'', '20', '', '30'};
% valence_label_1 = {'1', '2', '3', '4', '5'};
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};

% idxTheta = find(LF_steps >= 4 & LF_steps <= 8);
idxTheta = find(LF_steps >= 13 & LF_steps <= 30);
f = figure('units', 'centimeters', 'position', [5 5 3 2.5], 'Color', [1, 1, 1]); hold on;
h = cell(n_valence, 1);
comdlgrm_avg_valence = cell(n_valence, 1);
for j = 1:n_valence
    comdlgrm_avg_1 = mean(comdlgrm_valence{j, i_start_a, i_start_a}(:, idxTheta, :), 3);
    h{j} = plot(LF_steps, mean(comdlgrm_avg_1, 1), 'Color', color_valence(j, :), 'LineWidth', 0.5); hold on;
    ci = bootci(1000, @(x) mean(x), comdlgrm_avg_1);
    fill([LF_steps fliplr(LF_steps)], [ci(2,:) fliplr(ci(1,:))], color_valence(j, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    comdlgrm_avg_valence{j} = comdlgrm_avg_1;
end
xlim([13 30])
ylabel('MI');
xlabel('Beta frequency')
% ytickformat('%.1f'); 
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
    'YLim', [0.00002 0.0002], 'YTick', YTickRange, ...
    'XLim', [13 30], 'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01]);
% lgd = legend([h{1}, h{2}, h{3}, h{4}, h{5}], valence_label_1); legend('boxoff');
% lgd.FontSize = 3;
% lgd.ItemTokenSize = [6, 18];

outputFigure(fullfile([path dataName(1:4)]), [dataName '_stat']);

%--- save data for ART ANOVA in R
data1 = comdlgrm_avg_valence{1}; % 30 trials x 21 timepoints
data2 = comdlgrm_avg_valence{2}; % 120 trials x 21 timepoints
data3 = comdlgrm_avg_valence{3}; % 30 trials x 21 timepoints
data4 = comdlgrm_avg_valence{4}; % 30 trials x 21 timepoints
data5 = comdlgrm_avg_valence{5}; % 30 trials x 21 timepoints
time = LF_steps; % 1x 21 timepoints
save([dataName(1:4) '\' dataName '_PAC.mat'], "data1", "data2", "data3", "data4", "data5", "time");


%--- Trial-level permutation test, time resolved
offset = [0.001:0.01:0.1]; lindWidth = 1; color_all = repmat([1 1 1], 10, 1);
lab = yticklabels;
labV = zeros(size(lab, 1),1);
for iLab = 1:length(labV)
    labV(iLab) = str2num(lab(iLab, :));
end
time = LF_steps;
data_list = {comdlgrm_avg_valence{1}, comdlgrm_avg_valence{2}, comdlgrm_avg_valence{3}, comdlgrm_avg_valence{4}, comdlgrm_avg_valence{5}};  % cell array of 3 datasets
[pairs, sig_timepoints_all] = timeResolved_permutation_multiple_conditions(data_list, time, labV, lindWidth, offset, color_all, 0);


function [pairs, sig_timepoints_all] = timeResolved_permutation_multiple_conditions(data_list, time, labV, lindWidth, offset, color, visualMark)
n_perm = 1000;
alpha_cluster = 0.05;
pairs = nchoosek(1:length(data_list), 2);
n_pairs = size(pairs, 1);
n_time = size(data_list{1}, 2);
sig_timepoints_all = false(n_pairs, n_time);
all_cluster_pvals = {};  % store cluster p-values
all_cluster_labels = {};  % store cluster masks
perm_null_all = [];  % collect null cluster sums for FDR
for pair = 1:n_pairs
    idx1 = pairs(pair, 1);
    idx2 = pairs(pair, 2);
    d1 = data_list{idx1};
    d2 = data_list{idx2};
    n1 = size(d1, 1);
    n2 = size(d2, 1);
    % Compute observed t-values
    m1 = mean(d1, 1, 'omitnan');
    m2 = mean(d2, 1, 'omitnan');
    s1 = std(d1, 0, 1, 'omitnan');
    s2 = std(d2, 0, 1, 'omitnan');
    t_obs = (m1 - m2) ./ sqrt((s1.^2 / n1) + (s2.^2 / n2));
    % Cluster-forming threshold
    t_thresh = tinv(1 - alpha_cluster/2, min(n1,n2) - 1);
    cluster_mask = abs(t_obs) > t_thresh;
    [clust_labels, n_clust] = bwlabel(cluster_mask);
    cluster_sums_obs = zeros(1, n_clust);
    for i = 1:n_clust
        cluster_sums_obs(i) = sum(abs(t_obs(clust_labels == i)));
    end
    % Build permutation null
    combined = [d1; d2];
    perm_clust_max = zeros(n_perm, 1);
    for p = 1:n_perm
        perm_idx = randperm(n1 + n2);
        pd1 = combined(perm_idx(1:n1), :);
        pd2 = combined(perm_idx(n1+1:end), :);
        pm1 = mean(pd1, 1, 'omitnan');
        pm2 = mean(pd2, 1, 'omitnan');
        ps1 = std(pd1, 0, 1, 'omitnan');
        ps2 = std(pd2, 0, 1, 'omitnan');
        t_perm = (pm1 - pm2) ./ sqrt((ps1.^2 / n1) + (ps2.^2 / n2));
        mask = abs(t_perm) > t_thresh;
        [lbl, n_c] = bwlabel(mask);
        cluster_sums = zeros(1, n_c);
        for k = 1:n_c
            cluster_sums(k) = sum(abs(t_perm(lbl == k)));
        end
        perm_clust_max(p) = max([cluster_sums, 0]);
    end
    % Store for FDR correction
    perm_null_all = [perm_null_all; perm_clust_max];
    % Get cluster-level p-values
    pvals = zeros(1, n_clust);
    for i = 1:n_clust
        pvals(i) = mean(perm_clust_max >= cluster_sums_obs(i));
    end
    % Save
    all_cluster_pvals{pair} = pvals;
    all_cluster_labels{pair} = clust_labels;
end
% FDR correction over all cluster p-values
flat_pvals = cell2mat(all_cluster_pvals(:)');
fdr_pvals = mafdr(flat_pvals, 'BHFDR', true);
% Assign corrected p-values back
idx = 1;
for pair = 1:n_pairs
    this_labels = all_cluster_labels{pair};
    n_clust = max(this_labels);
    sig_timepoints = false(1, n_time);
    for c = 1:n_clust
        if fdr_pvals(idx) < 0.05
            sig_timepoints(this_labels == c) = true;
        end
        idx = idx + 1;
    end
    sigline = nan(1,numel(time));
    % idx = burstWindow_2(1):burstWindow_2(end);
    sigline(sig_timepoints) = labV(1) + offset(pair);
    if visualMark
        hold on
        plot(time,sigline, '-', 'LineWidth', lindWidth, 'Color', color(pair, :))
    end
    sig_timepoints_all(pair, :) = sig_timepoints;
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