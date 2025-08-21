clear; close all; clc;
dbstop if error;

path = 'F:\OneDrive\OneDrive - Yale University\Yajun\IAPS\Manuscript\Figures\';
%% Fig. 1c, Distribution of emotional rating
clear -path; close all;
load([path '\fig1\fig1.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
close all;
y = mean(n_rating, 2);
percent = y/sum(y)*100 
y_std = std(n_rating, 0, 2);
percent_std = y_std/sum(y)*100
sliceLabels = arrayfun(@(x) sprintf('%.1f%%', x), percent, 'UniformOutput', false);
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6], 'Color', [1, 1, 1]);
fontNumbleTick = 6; LineWidth = 0.25;
hPie = pie(percent, sliceLabels);
textHandles = findobj(hPie, 'Type', 'text');
patchHandles = findobj(hPie, 'Type', 'Patch');   
for i = 1:numel(textHandles)
    pos = get(textHandles(i), 'Position');   % [x, y, z]
    set(textHandles(i), 'Position', 0.5 * pos, 'FontSize', 6);  % Halfway towards the origin
    set(textHandles(i), 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment',   'middle');
    set(patchHandles(i), 'FaceColor', color_valence(i,:), 'LineWidth', 1.5, 'FaceAlpha', 0.6);
end 
set(patchHandles, 'LineWidth', 0.5);
% title('The proportions of valence rating');
lgd = legend(valence_label_1, 'Location','eastoutside', FontSize=5); legend('boxoff');
% lgd = legend([p1, p2, p3], {'5', '2 3 4', '1'}, 'box', 'off', 'fontsize', fontNumbleTick, 'units', 'centimeters');
lgd.FontSize = 6; 
lgd.ItemTokenSize = [10, 18];

figName = 'fig1c';
outputFigure(fullfile([path figName(1:4)]), figName);

%% Fig. 1d, Response time for emotional rating
clear -path; close all;
load([path '\fig1\fig1.mat' ])
valence_label_1 = {'VP', 'P', 'N', 'U', 'VU'};
figName = 'fig1d';
threshMin = 0.2; threshMax = 4.5;
f = figure('units', 'centimeters', 'position', [5, 5, 6, 5], 'Color', [1, 1, 1]);
fontNumbleTick = 6; LineWidth = 0.25;
hold on;
data_points = vertcat(RTs_valence{:});
idx = find(data_points >= threshMin & data_points <= threshMax);
data_points = data_points(idx); % remove outliers
lengths = cellfun(@length, RTs_valence);
labels = arrayfun(@(k) repmat({num2str(k)}, lengths(k), 1), 1:numel(lengths), 'UniformOutput', false);
labels = vertcat(labels{:});
labels = labels(idx);
vs = violinplot(data_points, labels, 'MedianMarkerSize', 15, 'MarkerSize', 0.25);
for i = 1:length(vs)
    vs(i).ShowData = true;
    vs(i).ShowMean = true;
    vs(i).ViolinColor = {color_valence(i, :)};
    vs(i).EdgeColor = [0.3, 0.3 , 0.3];
    vs(i).BoxColor = [0.3, 0.3 , 0.3];
    % vs(i).MedianMarkerSize = 0.5;
    vs(i).ViolinPlot.LineWidth = 0.25;
    vs(i).WhiskerPlot.LineWidth = 0.25;
end
ylim([0 5]);
ylabel('Response time (s)');
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', 1:2:5, 'YTickLabel', 1:2:5,'XTick', 1:5, 'XTickLabel', valence_label_1, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]); box off;

%--- save data for ART ANOVA in R
all_data = [];
all_conditions = [];
for i = 1:length(RTs_valence)
    n = numel(RTs_valence{i});
    all_data = [all_data; RTs_valence{i}];
    all_conditions = [all_conditions; repmat(i, n, 1)];
end
% Save to .mat file
save([figName(1:4) '\' figName '_RTs_valence.mat'], 'all_data', 'all_conditions');

%--- Trial-level permutation test
num_permutations = 1000;
pairs = nchoosek(1:5, 2);
nPairs = size(pairs, 1);
% nSamples = 30;
p_vals = nan(nPairs,1);
for i = 1:nPairs
    c1 = pairs(i,1);
    c2 = pairs(i,2);
    data1 = RTs_valence{c1};  % condition i
        data2 = RTs_valence{c2};  % condition j
        n1 = length(data1);
        n2 = length(data2);
        % Compute observed mean difference
        mean1 = mean(data1, 'omitnan');
        mean2 = mean(data2, 'omitnan');
        std1 = std(data1, 0, 'omitnan');
        std2 = std(data2, 0, 'omitnan');
        observed_t = (mean1 - mean2) / sqrt(std1^2/n1 + std2^2/n2);
        % Permutation
        combined_data = [data1; data2];
        perm_t = zeros(num_permutations, 1);
        for p = 1:num_permutations
            shuffled_idx = randperm(n1 + n2);
            perm_data1 = combined_data(shuffled_idx(1:n1));
            perm_data2 = combined_data(shuffled_idx(n1+1:end));
            m1 = mean(perm_data1, 'omitnan');
            m2 = mean(perm_data2, 'omitnan');
            s1 = std(perm_data1, 0, 'omitnan');
            s2 = std(perm_data2, 0, 'omitnan');
            perm_t(p) = (m1 - m2) / sqrt(s1^2/n1 + s2^2/n2);
        end
        % Two-tailed p-value
        p_val = sum(abs(perm_t) >= abs(observed_t)) / num_permutations;
       p_vals(i)  = p_val;
end
pvals_fdr = mafdr(p_vals, 'BHFDR', true);
sig_timepoints = find(pvals_fdr < 0.05);
disp('pairs of significance: ');
disp(num2str(pairs(sig_timepoints, :)));
for j = 1:length(sig_timepoints)
H = sigstar({valence_label_1(pairs(sig_timepoints(j),:))}, pvals_fdr(sig_timepoints(j)));
end

outputFigure(fullfile([path figName(1:4)]), figName);

%% Fig. 1e, Pupil diameter response for emotional rating
clear -path; close all;
load([path '\fig1\fig1e.mat' ])
f = figure('units', 'centimeters', 'position', [5, 5, 6, 5], 'Color', [1, 1, 1]);
fontNumbleTick = 6; LineWidth = 1;
% axes('units', 'centimeters', 'position', [1.5, 0.9, 4.25, 3]);
hold on;
box off
ci = bootci(1000, @(x) mean(x,'omitnan'), Pupil_valence_1_all);
fill([Pupil_presentation_time fliplr(Pupil_presentation_time)], [ci(2,:) fliplr(ci(1,:))], color_valence(5, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p1 = plot(Pupil_presentation_time, nanmean(Pupil_valence_1_all, 1), 'Color', color_valence(5, :), 'LineWidth',LineWidth, 'LineStyle','-');
% ci = bootci(1000, @(x) mean(x,'omitnan'), Pupil_valence_2_all);
% fill([Pupil_presentation_time fliplr(Pupil_presentation_time)], [ci(2,:) fliplr(ci(1,:))], [102 204 0]./255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% p2 = plot(Pupil_presentation_time, nanmean(Pupil_valence_2_all, 1), 'Color', [102 204 0]./255, 'LineWidth',LineWidth);
ci = bootci(1000, @(x) mean(x,'omitnan'), Pupil_valence_3_all);
fill([Pupil_presentation_time fliplr(Pupil_presentation_time)], [ci(2,:) fliplr(ci(1,:))], [102 204 0]./255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(Pupil_presentation_time, nanmean(Pupil_valence_3_all, 1), 'Color', [102 204 0]./255, 'LineWidth',LineWidth, 'LineStyle',':');
ci = bootci(1000, @(x) mean(x,'omitnan'), Pupil_valence_4_all);
fill([Pupil_presentation_time fliplr(Pupil_presentation_time)], [ci(2,:) fliplr(ci(1,:))], color_valence(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p3 = plot(Pupil_presentation_time, nanmean(Pupil_valence_4_all, 1), 'Color', color_valence(1, :), 'LineWidth',LineWidth, 'LineStyle','--');
lab = yticklabels;
labV = zeros(length(yticklabels),1);
for iLab = 1:length(yticklabels)
    labV(iLab) = str2num(lab{iLab});
end
legend('boxoff');
xlabel('Time (s)')
ylabel('Pupil diameter')
xline(0, '--', 'LineWidth', 0.5, 'Color', [0, 0, 0]);
% lgd = legend([p1, p2, p3], {'5', '2 3 4', '1'}, 'box', 'off', 'fontsize', fontNumbleTick, 'units', 'centimeters');
lgd = legend([p1, p2, p3], {'VU', 'P/N/U', 'VP'}, 'box', 'off', 'fontsize', fontNumbleTick, 'units', 'centimeters');
lgd.FontSize = 6; 
lgd.ItemTokenSize = [6, 18];
set(gca, 'FontSize', fontNumbleTick, 'LineWidth', 0.25, 'XColor', 'k', 'YColor', 'k', ...
    'YTick', [-1 0 1], 'YTickLabel', [-1 0 1],'XTick', -1:3, 'XTickLabel', -1:3, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]); box off;

figName = 'fig1e';
outputFigure(fullfile([path figName(1:4)]), figName);

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
