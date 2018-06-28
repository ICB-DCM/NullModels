function  [LoewePredictions, TallaridaUBPredictions, TallaridaLBPredictions, HandPredictions, BlissPredictions, HSAPredictions] = CorrPlots(CellLine)
%% Corrplot
%
% Plots a subplot matrix containing the histogram of the single
% predicitions, the histogram of the difference between the predictions and
% a scatterplot of the predictions of all null models.
%
% Jakob

            LoeweRelativeIndizes = [];
            LoeweIndizes = [];
            
            BlissIndizes = [];
            HandIndizes = [];
            HSAIndizes = [];
            
            %Prediction vectors
            
            LoewePredictions = [];
            BlissPredictions = [];
            HandPredictions = [];
            HSAPredictions = [];
            TallaridaLBPredictions = [];
            TallaridaUBPredictions = [];
                        
            CombiResponses = [];


    for j = 1:length(CellLine.Combinations)
            
        
            LoeweRelativeIndizes = [LoeweRelativeIndizes; CellLine.Combinations{j}.LoeweRelativeIndex];
            LoeweIndizes= [LoeweIndizes; CellLine.Combinations{j}.LoeweIndex];
            
            BlissIndizes= [BlissIndizes; CellLine.Combinations{j}.BlissIndex];
            HandIndizes = [HandIndizes; CellLine.Combinations{j}.HandIndex];
            HSAIndizes = [HSAIndizes; CellLine.Combinations{j}.HSAIndex];
            
            %Prediction vectors
            
            LoewePredictions = [LoewePredictions; CellLine.Combinations{j}.LoewePrediction'];
            BlissPredictions = [BlissPredictions; CellLine.Combinations{j}.BlissPrediction'];
            HandPredictions = [HandPredictions; CellLine.Combinations{j}.HandPrediction'];
            HSAPredictions = [HSAPredictions; CellLine.Combinations{j}.HSAPrediction'];
            TallaridaLBPredictions = [TallaridaLBPredictions; CellLine.Combinations{j}.TallaridaPrediction(:, 1)];
            TallaridaUBPredictions = [TallaridaUBPredictions ; CellLine.Combinations{j}.TallaridaPrediction(:, 2)];
                        
            CombiResponses = [CombiResponses; CellLine.Combinations{j}.Response];
            
    end
    
    %% Output: Correlations between different Predictors/Indizes!
    
    meanResponse = nanmean(CombiResponses, 2);
    
    disp('Correlation of the different Predictors (Loewe, TLB, TUB, Hand, Bliss, HSA)');


    A = corrcoef([LoewePredictions, TallaridaLBPredictions, TallaridaUBPredictions, HandPredictions,  BlissPredictions, HSAPredictions], 'Rows', 'complete');
    disp(min(LoewePredictions-HandPredictions));
    sum(LoewePredictions<HandPredictions)
    
    %disp('Correlation of the different Indizes (Loewe, relative Loewe, Bliss, Hand, HSA)')
 
%% Plot

% corrplot([meanResponse, LoewePredictions, BlissPredictions, HandPredictions, HSAPredictions], 'varnames', {'Measurements', 'Loewe', 'Bliss', 'Hand', 'HSA'});
% 
% % Change the title of the figure to the name of the cell line
% fig = gcf;
% set(fig,'Name',CellLine.Name,'NumberTitle','off');

%% Plot: 
%
% 



    Names = {"Loewe", "TallaridaUB", "TallaridaLB", "Hand", "Bliss", "HSA"};
    values = [LoewePredictions, TallaridaUBPredictions, TallaridaLBPredictions, HandPredictions, BlissPredictions, HSAPredictions];
    
    values = ones(size(values)) - values; % Transform form 1->0 to 0->1
    
    
    alpha = 0.015; % The transparency value.
    %alpha=1;
    
    gap = .05* .8/6; % size of gap between plots;
    
    fontsize = 8;
    
    
    %scatter(HSAPredictions, HandPredictions, 100*ones(size(HSAPredictions)), 'Marker', '.', 'MarkerEdgeAlpha', aplha)
    
    figure('position', [120 42, 650, 650])
    
    
    for i=1:6
        
        s{i, i} = subplot('Position', [0.1+(i-1)*0.8/6 + gap, 0.1+(6-i)*0.8/6 + gap, 0.8/6 - 2* gap, 0.8/6 - 2* gap]);
        histogram(values(:, i), 'Normalization', 'probability');
        xlim([0 1])
        ylim([0 .3])
        set(gca,'XTickLabel',[], 'YTickLabel',[], 'TickLength', [0 0]);
        set(gca, 'XColor', 'b', 'YColor', 'k')
        %scatter(1-HSAPredictions, 1-HandPredictions, 220*ones(size(HSAPredictions)), 'Marker', '.', 'MarkerEdgeAlpha', aplha)
    
    end
    
    for i = 1:6
        for j=i+1:6
                
            % below diagonal
            s{i, j} = subplot('Position', [0.1+(j-1)*0.8/6 + gap, 0.1+(6-i)*0.8/6 + gap, 0.8/6 - 2 * gap, 0.8/6- 2 * gap]); % Atention here the position is given as [left hight ... ....] = [col row ... ...]
            
            %scatter(values(:, j), values(:, i), 100*ones(size(HSAPredictions)), 'Marker', '.', 'MarkerEdgeAlpha', alpha);
            histogram(values(:, i) - values(:, j), [-.4:0.05:.4], 'Normalization', 'probability', 'FaceColor', 'r')
            hold on
            line([0 0], [0 1], 'LineStyle','--', 'Color', 'k', 'LineWidth', 2);
            xlim([-.4 .4])
            set(gca, 'XAxisLocation','top', 'YAxisLocation','right', 'XTickLabel',[], 'YTickLabel',[], 'TickLength', [0.05 0.05])
            set(gca, 'XColor', 'r', 'YColor', 'r') % Change color Values of the axes
            
            %xlim([0 1])
            ylim([0 1])
            %set(gca,'XTickLabel',[], 'YTickLabel',[], 'TickLength', [0.05 0.05]);
            box on
            
            
            
            %above diagonals
            s{j, i} = subplot('Position', [0.1+(i-1)*0.8/6 + gap, 0.1+(6-j)*0.8/6 + gap, 0.8/6 - 2 * gap, 0.8/6 - 2 * gap]);
            scatter(values(:, i), values(:, j), 100*ones(size(HSAPredictions)), 'Marker', '.', 'MarkerEdgeAlpha', alpha);
            %xlim([0 1])
            %ylim([0 1])
            set(gca,'XTickLabel',[], 'YTickLabel',[], 'TickLength', [0.05 0.05]);
            set(gca, 'XColor', 'b', 'YColor', 'b')
            box on
            
        end
    end
    
    % Fixing the upper left Axis-Label
    set(gcf, 'currentaxes', s{1, 1})
    set(gca, 'XAxisLocation','top') 
    
    %asign temporal xticks to get the position of the xlabel
    xticks([0 0.5 1]), xticklabels([0 0.5 1]);
    xlabel(Names{1})
    
    yticks([0 0.15 .3]), yticklabels([0 0.5 1]);
    ylabel(Names{1})
    
    LoeweXLabel = get(gca, 'XLabel');
    XLabelPosition = LoeweXLabel.Position;

    LoeweYLabel = get(gca, 'YLabel');
    YLabelPosition = LoeweYLabel.Position;
    
    xticklabels([]);
    yticklabels([]);
    % Asign the new Position after the next for loop
    
    
    for i = 1:6
        
        set(gcf, 'currentaxes', s{1, i}) % the upper row
        xlabel(Names{i}, 'Interpreter', 'tex', 'Fontsize', fontsize, 'Color','k');
        
        set(gcf, 'currentaxes', s{6, i}) % the lowest row
        xticks([0 0.5 1]);
        xticklabels([0 0.5 1]);
        
        set(gcf, 'currentaxes', s{i, 1}) % the left col
        yticks([0 0.5 1]);
        yticklabels([0 0.5 1]);
        ylabel(Names{i}, 'Interpreter', 'tex', 'Fontsize', fontsize, 'Color','k');
        
    end
    
    
    set(gcf, 'currentaxes', s{1, 1})
    yticklabels([]);
    set(LoeweXLabel, 'Position', XLabelPosition);
    set(LoeweYLabel, 'Position', YLabelPosition);
	
    %% Now the x and yticks of the upper and right row have to be adjusted to the "distribution of the difference"-Plot
    
    for i = 1:5
        
            set(gcf, 'currentaxes', s{1, i+1})
            set(gca, 'xtick', [-.3, 0, .3], 'xticklabel', [-.3, 0, .3])
            
            set(gcf, 'currentaxes', s{i, 6})
            set(gca, 'ytick', [0, 0.5, 1], 'yticklabel', [0, 0.5, 1])
            
    end
    
    set(gcf, 'currentaxes', s{1, 1})
    yticks([0 0.15 .3]) 
    yticklabels([0 0.5 1]);
    
    for i = 1:6
        for j = i+1:6
            set(gcf, 'currentaxes', s{j, i});
            text(.07, .85, num2str(A(i, j), '%4.3f'), 'Interpreter', 'tex', 'Fontsize', fontsize);
            %above
            set(gcf, 'currentaxes', s{i, j});
            text(-.35, .85, num2str(A(i, j), '%4.3f'), 'Interpreter', 'tex', 'Fontsize', fontsize);
        end
    end
    
    % adjust yticks of Plot(1, 1)
    set(gcf, 'currentaxes', s{1, 1})
    yticklabels([]);
    
    end

%% Outtakes:

%% Indizes Histogramms
    
%     figure()
%     title('Histogram of Combination Indizes')
%     subplot(2, 2, 1)
%     hist(LoeweIndizes)
%     title('Loewe Indizes')
% 
%     subplot(2, 2, 2)
%     hist(BlissIndizes)
%     title('Bliss Indizes')
% 
%     subplot(2, 2, 3)
%     hist(HandIndizes)
%     title('Hand Indizes')
% 
%     subplot(2, 2, 4)
%     hist(HSAIndizes)
%     title('HSA Indizes')
%     
%     figure()
%     subplot(2, 2, 1)
%     plot(LoewePredictions, BlissPredictions, 'x')
%     xlabel('Loewe Predictions')
%     ylabel('Bliss Predictions')
% 
%     subplot(2, 2, 2)
%     plot(LoewePredictions, HandPredictions, 'x')
%     xlabel('Loewe Predictions')
%     ylabel('Hand Predictions')
% 
%     subplot(2, 2, 3)
%     plot(LoewePredictions, HSAPredictions, 'x')
%     xlabel('Loewe Predictions')
%     ylabel('HSA Predictions')
% 
%     subplot(2, 2, 4)
%     spy()

