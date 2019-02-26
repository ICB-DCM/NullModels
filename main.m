%% 
%
% This script reads in the data from O'Neil et al.


singleDrugData = 'FILL_IN_YOUR_DIRECTORY_HERE\Data\156849_1_supp_0_w2lh45.xlsx';
% change to file path of the single drug data of O'Neil et al. 2016

CombinationData = 'FILL_IN_YOUR_DIRECTORY_HERE\Data\156849_1_supp_1_w2lrww.xls';
% change to file path of the combination data of O'Neil et. al 2016

%% Read in data
tic()
  D = Data(singleDrugData, CombinationData);
  disp('Time for reading in the Data + Fitting the Hill-Curves')    
toc()


%% evaluate the synergy and print out computation time

times = zeros(5, 1); % Stores the computation times for the Synergy-Evaluation

 for i = 1:length(D.CellLines)
    
    %% Loewe
    
	tic()
        D.CellLines{i}.evaluateSynergyLoewe();
    times(1) = times(1) + toc();
    
    %% Bliss
    
    tic()
        D.CellLines{i}.evaluateSynergyBliss();
    times(2) = times(2) + toc();
    
    %% Hand
    
    tic()
        D.CellLines{i}.evaluateSynergyHand();
    times(3) = times(3) + toc();
    
    %% HSA
    tic()
        D.CellLines{i}.evaluateSynergyHSA();
    times(4) = times(4) + toc();

    %% Tallarida
    
    tic()
        D.CellLines{i}.evaluateSynergyTallarida();
	times(5) = times(5) + toc();
    
    disp(strcat(num2str(i), ' cell lines computed'))
    disp('computation times:')
    disp(times);
    %keyboard;

 end

disp('End of Data Evaluation')


%% Access the reliability of the Fits of the Single-Drug-Data:
%
% Check for each Drug fit at each cell line the residual, and if the
% predictions are valid (= $\in [0, 1]$). Additionally plot the Parameters
% of the drug fits

residuals = [];
noMeasurements = [];

validResponses = [];

ws = [];
ds = [];
ds_normalized = []; % = d / max_Dose
ns = [];


for i = 1:length(D.CellLines)
    for j = 1:length(D.CellLines{i}.Drugs)
        
        residuals = [residuals; D.CellLines{i}.Drugs{j}.FittingResidual];
        
        ws = [ws; D.CellLines{i}.Drugs{j}.Parameters(1)];
        ds = [ds; D.CellLines{i}.Drugs{j}.Parameters(2)];
        ds_normalized = [ds_normalized; D.CellLines{i}.Drugs{j}.Parameters(2)/max(D.CellLines{i}.Drugs{j}.Conc) ];
        ns = [ns; D.CellLines{i}.Drugs{j}.Parameters(3)];
        
        noMeasurements = [noMeasurements; length(D.CellLines{i}.Drugs{j}.Response)];
        validResponses = [validResponses; ~ (sum(D.CellLines{i}.Drugs{j}.Response <0) && sum(D.CellLines{i}.Drugs{j}.Response >1))];
    end
end


%% Plots
%
% Filter out the Drug-Fits with w_max==0, since those belong to the
% constant zero solution choosen by the BIC

figure()

% w_s
subplot(2, 2, 1)
histogram(ws(ws ~=0), 20, 'Normalization','probability')
title('w_{max}')
xlabel(strcat(num2str(sum(ws==0)), ' zero-response curves '))

% d_s
subplot(2, 2, 2)
histogram(ds_normalized(ws ~=0), 20, 'Normalization','probability')
title('d_{1/2} / maxDose')
xlim([0, 1])

% w_s
subplot(2, 2, 3)
histogram(ns(ws ~=0), 20, 'Normalization','probability')
title('Hill coeffitient \in [1, 20]')
xlim([1, 20])
xticks([1 2:2:20])

%avg residual
subplot(2, 2, 4)
histogram(residuals./noMeasurements , 20, 'Normalization','probability')
title('Residual/noMeasurement')

disp('Number of zero response fits:')
disp(sum(ws==0))


%% Reproduce the plots from the publication

IsobolePlot(D, 2, 3)

CorrPlots(D.CellLines{3});

VolumeMetricConceptPlot

VolPlot(D);