n_valence = 5;
responsive_marks = {'enhanced', 'suppressed'};
for respNeuronIdx = 1:n_valence

    for iRepMode = 1:length(responsive_mode)
        switch responsive_mode{iRepMode}
            case 'excitatory'
                ind_tmp = indices_valence.excitatory{respNeuronIdx};
            case 'inhibitory'
                ind_tmp = indices_valence.inhibitory{respNeuronIdx};
        end

        
        %% GLM
        FR_valence_all_tmp = []; label_val = [];
        for i_val = 1:n_valence
            FR_tmp = cell2mat(FR_valence_all{i_val}(ind_tmp)');
            label_tmp = repmat(i_val, size(FR_tmp, 1), 1);
            label_val = cat(1, label_val, label_tmp);
            FR_valence_all_tmp = cat(1, FR_valence_all_tmp, FR_tmp);
        end
        label_val = categorical(label_val);
        
        nanIdx = find(any(isnan(FR_valence_all_tmp), 2));
        FR_valence_all_tmp(nanIdx, :) = [];
        label_val(nanIdx) = [];

        [~, idx_1] = min(abs(trialTime - (0.3)));
        [~, idx_2] = min(abs(trialTime - (1.3)));
        tw_task = idx_1:idx_2;
        [~, idx_1] = min(abs(trialTime - (-1.0)));
        [~, idx_2] = min(abs(trialTime - (0)));
        tw_bs = idx_1:idx_2;
        baselineType = 'normchange'; % 'absolute', 'relative', 'relchange', 'normchange', 'db', 'zscore', 'no'
        tw_norm = find(trialTime >= -1 & trialTime <= 3);
        FR_valence_all_tmp_1 = firingRate_normalization(FR_valence_all_tmp, baselineType, 1:size(FR_valence_all_tmp, 2));
        FR_task = mean(FR_valence_all_tmp_1(:, tw_task), 2);
        FR_bl = mean(FR_valence_all_tmp_1(:, tw_bs), 2);
       
        if respNeuronIdx == 1
            label_val1 = label_val;
        elseif respNeuronIdx == 2
            label_val1 = reordercats(label_val, {'2','1','3','4','5'});
        elseif respNeuronIdx == 3
            label_val1 = reordercats(label_val, {'3','1','2','4','5'});
        elseif respNeuronIdx == 4
            label_val1 = reordercats(label_val, {'4','1','2','3','5'});
        elseif respNeuronIdx == 5
            label_val1 = reordercats(label_val, {'5','1','2','3','4'});
        end

        data_tbl = table(FR_task, FR_bl, label_val1, 'VariableNames', {'FiringRate_task', 'BaselineFR', 'Valence'});
        % mdl3 = fitglm(data_tbl, 'FiringRate_task ~ BaselineFR + Valence', 'Distribution', 'poisson', 'Link', 'log');
        mdl3 = fitglm(data_tbl, 'FiringRate_task ~ BaselineFR + Valence', 'Distribution', 'normal');

        residual = mdl3.Residuals(:, 1);
        MSE = table2array(mean(residual.^2));
        R_squared = mdl3.Rsquared.Ordinary;
        VIF = 1/(1 - R_squared^2);
        beta_estimate  = mdl3.Coefficients.Estimate';

        residuals = mdl3.Residuals.Raw;  % or use .Standardized for normalized residuals
        data_tbl.Residuals = residuals;
        grouped_resid = groupsummary(data_tbl, "Valence", "mean", "Residuals");
        disp(grouped_resid);


        p_vals = mdl3.Coefficients.pValue(contains(mdl3.CoefficientNames, 'Valence'));
        pvals_fdr = mafdr(p_vals, 'BHFDR', true);
        disp( [dataName ', p values for Pref_' num2str(valence_number(respNeuronIdx)) ', ' responsive_marks{iRepMode} ', compared valence: ' num2str(setdiff(1:5, respNeuronIdx))])
        disp(num2str(pvals_fdr));

        %--- plot the beta coefficients
        coefs = mdl3.Coefficients
        b = coefs.Estimate(2:end);
        ci = coefCI(mdl3); ci = ci(2:end, :);
        close all;
        f = figure('units', 'centimeters', 'position', [5 5 4 3], 'Color', [1, 1, 1]); hold on;
        box off
        xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);
        % xlabel('Time (s)', 'FontSize', 6);
        ylabel('Coefficient Estimate', 'FontSize', 6)
        set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
            'XLim', [0 5], 'XTick', 1:5, 'XTickLabel', {'FR\_B', valence_label_1{setdiff(1:5, respNeuronIdx)}},...
            'YLim', [-0.2 0.5], 'YTick', -0.2:0.2:0.4, 'YTickLabel', -0.2:0.2:0.4, ...
            'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
        errorbar(1:length(b), b, b - ci(:,1), ci(:,2) - b, 'o', 'Color', [0, 0, 0], 'MarkerSize', 4, 'LineWidth', 0.5);

        outputFigure(fullfile([path dataName(1:4)]), [dataName '_Pref' num2str(valence_number(respNeuronIdx)) '_' responsive_marks{iRepMode} '_GLM_Coefficient']);


        %--- plot the predicted FRs
        valence_levels = categories(data_tbl.Valence);
        mean_baseline = mean(data_tbl.BaselineFR);
        % Create new data table for prediction
        new_data = table(repmat(mean_baseline, numel(valence_levels), 1), ...
            categorical(valence_levels, valence_levels), ...
            'VariableNames', {'BaselineFR', 'Valence'});
        % Predict fitted firing rate
        [yhat, yci] = predict(mdl3, new_data);
        close all;
        f = figure('units', 'centimeters', 'position', [5 5 4 3], 'Color', [1, 1, 1]); hold on;
        box off
        xline(0, '--', 'LineWidth', 0.25, 'Color', [0, 0, 0]);
        % xlabel('Time (s)', 'FontSize', 6);
        ylabel('Predicted FRs', 'FontSize', 6)
        set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
            'XLim', [0 5], 'XTick', 1:5, 'XTickLabel', valence_label_1(str2double(valence_levels)),...
            'YLim', [0.2 0.4], 'YTick', 0.2:0.1:0.4, 'YTickLabel', 0.2:0.1:0.4, ...
            'tickdir', 'out', 'ticklength', [0.01, 0.01]); box off;
        errorbar(1:numel(valence_levels), yhat, yhat - yci(:,1), yci(:,2) - yhat, 'o', 'Color', [0, 0, 0], 'MarkerSize', 4, 'LineWidth', 0.5);

        outputFigure(fullfile([path dataName(1:4)]), [dataName '_Pref' num2str(valence_number(respNeuronIdx)) '_' responsive_marks{iRepMode} '_GLM_PredictedFR']);


    end

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