%% VolumeMetricConceptPlot
%
% In this script a visualisation of the Volume-Score for synergy detection
% is given. Herefore in one figure the (averaged) measured response as well
% as the predicted response curve is plotted. In a second the volume of
% intersection is plotted.
%
% Jakob


%% Color Scheme:
%
% Red/ Antagonistic:    [166/256, 94/256, 185/256]
% Yellow/ Neutral:          [240/256, 240/256, 50/256]
% Green/ Synergistic:   [59/256, 195/256, 123/256]



i=1; % Cell line number: Change here
j=13; % Combination number: Change here 

C = D.CellLines{i}.Combinations{j}; % C is the Combination, that is used

CellLineName = D.CellLines{i}.Name;


% Get the axis values 
ConcA = log(C.ConcA);
ConcB = log(C.ConcB);

%Shift
ConcA = (ConcA-min(ConcA)) / (max(ConcA)-min(ConcA)); 
ConcB = (ConcB-min(ConcB)) / (max(ConcB)-min(ConcB));


% Plot the measurement data
figure()
for i = 1:4
    %scatter3(ConcA, ConcB, C.Response(:, i), 'filled', 'MarkerFaceColor', 'b')
    scatter3(ConcA, ConcB, 1-C.Response(:, i), 'filled', 'MarkerFaceColor', 'b') % 1- Response to get the Maximal response to 1.
    hold on
end

%scatter3(ConcA, ConcB, nanmean(C.Response, 2), 'filled', 'MarkerFaceColor',[0 .75 .75])


AxisA = unique(ConcA, 'sorted');
AxisB = unique(ConcB, 'sorted');

% xy plane axes
[X, Y] = meshgrid(AxisA, AxisB);
Z1 = zeros(size(X));
Z2 = zeros(size(X));


% Color of the surfaces
colorRed = [166/256, 94/256, 185/256];
colorYellow = [240/256, 240/256, 50/256];
colorGreen = [59/256, 195/256, 123/256];

C1(:, :, 1) = colorYellow(1) * ones(size(X));
C1(:, :, 2) = colorYellow(2) * ones(size(X));
C1(:, :, 3) = colorYellow(3) * ones(size(X));

C2(:, :, 1) = colorGreen(1) *ones(size(X));
C2(:, :, 2) = colorGreen(2)*ones(size(X));
C2(:, :, 3) = colorGreen(3) *ones(size(X));

C3(:, :, 1) = colorRed(1) *ones(size(X));
C3(:, :, 2) = colorRed(2)*ones(size(X));
C3(:, :, 3) = colorRed(3) *ones(size(X));


%% Plot mead response  surface

for i = 1:4
    for j=1:4
       
        %Z1(j, i) = nanmean(C.Response(intersect(find(ConcA == AxisA(i)), find(ConcB == AxisB(j))), :));
        Z1(j, i) = 1-nanmean(C.Response(intersect(find(ConcA == AxisA(i)), find(ConcB == AxisB(j))), :));
        
    end
end

surf(X, Y, Z1, C1, 'FaceAlpha', 0.6);
hold on

%% Plot Hand Prediction surface

for i = 1:4
    for j=1:4
       
        %Z2(j, i) = C.HandPrediction(intersect(find(ConcA == AxisA(i)), find(ConcB == AxisB(j))));
        Z2(j, i) = 1-C.HandPrediction(intersect(find(ConcA == AxisA(i)), find(ConcB == AxisB(j))));
        
    end
end

surf(X, Y, Z2, C2, 'FaceAlpha', 0.6);


% plot formatting

xlabel(strcat(C.DrugA.Name, ' dose'))
ylabel(strcat(C.DrugB.Name, ' dose'))
title(strcat('Cell line: ', CellLineName))

zlim([-0.1 0.8])
axis square
%view(125, 15)
view(-55, 15) % Angle on the plot

disp('Dots: Measurements, Red: (avg) Response Surface, Blue: Predicted Response Surface')


% Plot the polygone between the two surfaces in an additional plot

figure()
surf(X, Y, Z1, C1, 'FaceAlpha', 0.8);
hold on
surf(X, Y, Z2, C1, 'FaceAlpha', 0.8);

% Faces with constant x axis
hold on
surf(zeros(4, 2), Y(:, [1, 1]), [Z1(:, 1), Z2(:, 1)], C1(:, [1 2], :), 'FaceAlpha', 0.8)
hold on
surf(ones(4, 2), Y(:, [4, 4]), [Z1(:, 4), Z2(:, 4)], C1(:, [1 2], :), 'FaceAlpha', 0.8)
hold on

% Faces with constant y axis
hold on
surf(X([1, 1], :)', zeros(4, 2), [Z1(1, :)', Z2(1, :)'], C1(:, [1 2], :), 'FaceAlpha', 0.8)
hold on
surf(X([4, 4], :)', ones(4, 2), [Z1(4, :)', Z2(4, :)'], C1(:, [1 2], :), 'FaceAlpha', 0.8)
hold on


% Set all the axis limits and views to the same values as in the response
% surface plot.

zlim([-0.1 0.8])
axis square
%view(125, 15)
view(-55, 15)
set(gca, 'Visible', 'off') % make the axes invisible
