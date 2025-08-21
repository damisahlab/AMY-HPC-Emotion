%% FIg. S3, FRs of pyramidal and interneurons during beta bursts
figName_all = {'figS3a1_HPC_pyr', 'figS3a2_HPC_pyr', 'figS3a3_HPC_int', 'figS3a4_HPC_int',...
    'figS3b1_AMY_pyr', 'figS3b2_AMY_pyr', 'figS3b3_AMY_int', 'figS3b4_AMY_int',...
    'figS3c1_HPC-AMY_pyr', 'figS3c2_HPC-AMY_pyr', 'figS3c3_HPC-AMY_int', 'figS3c4_HPC-AMY_int',...
    'figS3d1_AMY-HPC_pyr', 'figS3d2_AMY-HPC_pyr', 'figS3d3_AMY-HPC_int', 'figS3d4_AMY-HPC_int'};
for i_fig = 1:length(figName_all)
    figName = figName_all{i_fig};
    clear -path; close all;
    fontNumbleLable = 6; fontNumbleTick = 6;
    LineWidth = 0.25;
    YTickRange = [0 0.5 1]; YTickLabel = [0 0.5 1]; XTickRange = [-0.2 0 0.2]; XTickLabel = [-0.2 0 0.2];
    fig = openfig([path 'figS3\' figName '.fig']); box off;
    ax = findall(fig, 'type', 'axes');
    for i = 1:length(ax)
        lines = findall(ax(i), 'type', 'line');
        for j = 1:length(lines)
            % lines(j).LineWidth = 0.5;
            if contains(lines(j).DisplayName, 'data') || isempty(lines(j).DisplayName)
                lines(j).LineWidth = 1.0;
                if strcmp(figName, 'fig4b1')
                    y = lines(i).YData;
                    lines(j).YData = y - 0.045;
                elseif strcmp(figName, 'fig4e1')
                    y = lines(i).YData;
                    lines(j).YData = y - 0.08;
                end
            else
                lines(j).LineWidth = 0.5;
            end
        end
        ax(i).Position = [0.18 0.21 0.7 0.8];
    end
    set(fig, 'Units', 'centimeters', 'Position', [4 4 4 3], 'Color', [1 1 1]);
    % axes('units', 'centimeters', 'position', [0.5, .9, 1.5, 2]); hold on;
    set(gca, 'FontSize', fontNumbleTick, 'LineWidth', LineWidth, 'XColor', 'k', 'YColor', 'k', ...
        'YTick', [0 0.1 0.2], 'YTickLabel', [0 0.1 0.2], 'XLim', [-0.2 0.2], 'XTick', XTickRange, 'XTickLabel', XTickLabel, ...
        'tickdir', 'out', 'ticklength', [0.01, 0.01]);
    title(gca, '');
    ylabel('Z-scored FR')
    xlabel('Time (s)')
    lgd = legend('Location','northeast'); legend('boxoff');
    lgd.FontSize = 6;
    % lgd.Position(1) = lgd.Position(1) + 0.25;
    lgd.ItemTokenSize = [6, 18];

    outputFigure(fullfile([path figName(1:5)]), figName);
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