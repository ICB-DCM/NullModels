function IsobolePlot(D, i, j)
%% IsobolePlot(D, i, j)
%
%   Makes a plot of the two dose Response Curves of the Drugs in
%   D.CellLines{i}.Combinations{j} and the Isoboles of the Corresponding
%   Combination.
%
%
%% List of Variables
%
%   C                                        The Combination
%
%   aDoseAxisLimits :           [1x2] Limits of the Axes corresp to DrugA
%   bDoseAxisLimits:                          -"-
%   ConcentrationA:               The Concentration axis values of the
%                                              DrugA-Dose-Resp-Plot
%   ResponsesA:                    The corresp Responses
%   ConcentrationB:                ...
%   ResponsesB:                    ...
%
%   DrugAPlot:                         The Subplot of DrugA
%   DrugBPlot:                         The Subplot of DrugB
%   IsoboleSubplot:                 The Subplot for the Isobboles
%
%   IsoboleResponse:            The Response of the Isobole
%   


% Parameters of the Plot

    %space = 'log';
    space = 'lin';

    figure()
    %set(gcf, 'Position', [100, 5, 695, 695]) % set the figure to be rectangular
    set(gcf, 'Position', [100, 5, 600, 600])
    
    C = D.CellLines{i}.Combinations{j};

    %aDoseAxisLimits = [0.1*min(C.DrugA.Conc), 10*max(C.DrugA.Conc)]; % The Limits of the Drug A Dose-Axis
    %bDoseAxisLimits = [0.1*min(C.DrugB.Conc), 10*max(C.DrugB.Conc)]; % The Limits of the Drug B Dose-Axis
    if strcmp(space, 'log')
        aDoseAxisLimits = [min(C.DrugA.Conc), 2*max(C.DrugA.Conc)]; % The Limits of the Drug A Dose-Axis
        bDoseAxisLimits = [min(C.DrugB.Conc), 2*max(C.DrugB.Conc)]; % The Limits of the Drug B Dose-Axis
    else
        aDoseAxisLimits = [0, max(C.DrugA.Conc)/3.7]; % The Limits of the Drug A Dose-Axis
        bDoseAxisLimits = [0, max(C.DrugB.Conc)/9.2]; % The Limits of the Drug B Dose-Axis
    end
    
    % If we want the same axis limits for both drugs
    %aDoseAxisLimits = [min(aDoseAxisLimits(1), bDoseAxisLimits(1)), max(aDoseAxisLimits(2), bDoseAxisLimits(2))];
    %bDoseAxisLimits = aDoseAxisLimits;
    
    % The Response-Line for the Isobole;
    
    if C.DrugA.Parameters(1)>0.5 && C.DrugB.Parameters(1)>0.5 % If 0.5 is in the dose range take that for the Isobole...
        IsoboleResponse = 0.5;
    else
        IsoboleResponse = 1- 1/2* min(C.DrugA.Parameters(1), C.DrugB.Parameters(1));
    end


%% Plot Drug A

    %Evaluate the Drug for plotting
    if strcmp(space, 'log')
        ConcentrationsA = log(aDoseAxisLimits(1)): (log(aDoseAxisLimits(2))-log(aDoseAxisLimits(1)))/100: log(aDoseAxisLimits(2)); % equidistant in Log-Space
        ConcentrationsA = exp(ConcentrationsA);
    else
        ConcentrationsA = aDoseAxisLimits(1) : (aDoseAxisLimits(2)-aDoseAxisLimits(1))/100: aDoseAxisLimits(2);
    end
    
        ResponsesA = zeros(size(ConcentrationsA));
    
    for i = 1:numel(ResponsesA)
        ResponsesA(i) = C.DrugA.evaluateDrug(ConcentrationsA(i));
    end

    % Plot
    DrugAPlot = subplot('Position', [0.1 0.3 0.18 0.65]);
    
    if strcmp(space, 'log')
    
        semilogy(1-ResponsesA, ConcentrationsA, 'LineWidth',2);
        hold on
        semilogy([0, 1.1], C.DrugA.invertDrug(IsoboleResponse)*[1,1], 'r --')
        hold on 
        semilogy(1-IsoboleResponse, C.DrugA.invertDrug(IsoboleResponse), 'r o')
        hold on 
        semilogy((1-IsoboleResponse)*[1, 1], aDoseAxisLimits, 'r--')
    
    else
        
        plot(1-ResponsesA, ConcentrationsA, 'LineWidth',2);
        hold on
        plot([0, 1.1], C.DrugA.invertDrug(IsoboleResponse)*[1,1], 'r --')
        hold on 
        plot(1-IsoboleResponse, C.DrugA.invertDrug(IsoboleResponse), 'r o')
        hold on 
        plot((1-IsoboleResponse)*[1, 1], aDoseAxisLimits, 'r--')

    end
        
    xlim([0 1.1])
    ylim(aDoseAxisLimits)
    ylabel(strcat(C.DrugA.Name, ' dose'))
    
    %xlabel('Response')
    
    set(gca, 'xdir', 'reverse')    

%% Plot Drug B
    
    if strcmp(space, 'log') 
        ConcentrationsB = log(bDoseAxisLimits(1)): (log(bDoseAxisLimits(2))-log(bDoseAxisLimits(1)))/100: log(bDoseAxisLimits(2));
        ConcentrationsB = exp(ConcentrationsB);
    else
        ConcentrationsB = bDoseAxisLimits(1): (bDoseAxisLimits(2)-bDoseAxisLimits(1))/100:bDoseAxisLimits(2);
    end
    
    ResponsesB = zeros(size(ConcentrationsB));

    for i = 1:numel(ResponsesB)
        ResponsesB(i) = C.DrugB.evaluateDrug(ConcentrationsB(i));
    end

    %DrugBPlot = subplot(2, 2, 4);
    drugBPlot = subplot('Position', [0.3 0.1 0.65 0.18]);
    
    if strcmp(space, 'log')
    
        semilogx(ConcentrationsB, 1-ResponsesB,'LineWidth',2);
        hold on
        semilogx(C.DrugB.invertDrug(IsoboleResponse)*[1,1], [0, 1.1], 'r --')
        hold on 
        semilogx(C.DrugB.invertDrug(IsoboleResponse), 1-IsoboleResponse, 'r o')
        hold on 
        semilogy(bDoseAxisLimits, (1-IsoboleResponse)*[1, 1], 'r--')
    
    else
        
        plot(ConcentrationsB, 1-ResponsesB,'LineWidth',2);
        hold on
        plot(C.DrugB.invertDrug(IsoboleResponse)*[1,1], [0, 1.1], 'r --')
        hold on 
        plot(C.DrugB.invertDrug(IsoboleResponse), 1-IsoboleResponse, 'r o')
        hold on 
        plot(bDoseAxisLimits, (1-IsoboleResponse)*[1, 1], 'r--')    

    end

    xlim(bDoseAxisLimits)
    ylim([0 1.1])
    xlabel(strcat(C.DrugB.Name, ' dose'))
    
    %h = ylabel('Response');
    %set(h,'Rotation',45);
    
    %% Insert Response- x/ylabel:
    
    subplot('Position', [0.1 0.1 0.15 0.15]);
    text(.4, .4, 'Response', 'Rotation', 45, 'FontSize', 11)
    axis off
    
%% Isoboles

    %IsoboleSubplot = subplot(2, 2, 2);
    IsoboleSubplot = subplot('Position', [0.3 0.3 0.65 0.65]);
    axis square


%% Loewe
    lambda = [0:1e-6:1];
    if strcmp(space, 'log')
        loglog( C.DrugB.invertDrug(IsoboleResponse) *lambda , C.DrugA.invertDrug(IsoboleResponse) * (1-lambda))
    else
        plot( C.DrugB.invertDrug(IsoboleResponse) *lambda , C.DrugA.invertDrug(IsoboleResponse) * (1-lambda))
    end


    %linkaxes([DrugAPlot, IsoboleSubplot], 'y')
    %linkaxes([DrugBPlot, IsoboleSubplot], 'x')
    xlim(bDoseAxisLimits)
    ylim(aDoseAxisLimits)

    set(IsoboleSubplot, 'xticklabels', [])
    set(IsoboleSubplot, 'yticklabels', [])

    %% Tallarida

    lambda = [1e-6, 1e-5, 1e-4, 1e-3, 0.01:0.01:0.99, 1-1e-3, 1-1e-4, 1-1e-5, 1-1e-6];
    %lambda =[1e-8, 1e-6, 1e-5, 0:1e-3:1, 1-1e-5, 1-1e-6, 1-1e-7, 1-1e-8, 1e-10];
    isoboleA = zeros(numel(lambda), 2);
    isoboleB = zeros(numel(lambda), 2);
    
    
    for i = 1:numel(lambda)

        aMin = lambda(i) * aDoseAxisLimits(1);
        aMax = 2*lambda(i) * aDoseAxisLimits(2);
        bMin =  (1-lambda(i)) * bDoseAxisLimits(1);
        bMax =2*(1-lambda(i)) * bDoseAxisLimits(2);
        
        [isoboleA(i, 1), isoboleB(i, 1)] = bisection(IsoboleResponse, @(a, b)C.tallaridaAB(a, b), aMin, bMin, aMax, bMax);
        [isoboleA(i, 2), isoboleB(i, 2)] = bisection(IsoboleResponse, @(a, b)C.tallaridaBA(a, b), aMin, bMin, aMax, bMax);

    end
    if strcmp(space, 'log')
        hold on
        loglog(max(isoboleB, [], 2), max(isoboleA, [], 2));
        hold on
        loglog(min(isoboleB, [], 2), min(isoboleA, [], 2));
    else
        hold on
        plot(max(isoboleB, [], 2), max(isoboleA, [], 2));
        hold on
        plot(min(isoboleB, [], 2), min(isoboleA, [], 2));
    end    
    
    
    
%% Hand
    
    isoboleA = zeros(numel(lambda), 1);
    isoboleB = zeros(numel(lambda), 1);
    
    
    for i = 1:numel(lambda)

        %aMin = lambda(i) * aDoseAxisLimits(1);
        aMin = 0;
        aMax = lambda(i) * aDoseAxisLimits(2);
        %bMin =  (1-lambda(i)) * bDoseAxisLimits(1);
        bMin = 0;
        bMax =(1-lambda(i)) * bDoseAxisLimits(2);
        
        if lambda(i) == 0
            
                isoboleA(i) = 0; 
                isoboleB(i) = IsoboleResponse;
        
        elseif lambda(i) ==1
            
                isoboleA(i) = IsoboleResponse; 
                isoboleB(i) = 0;
        
        else
            [isoboleA(i), isoboleB(i)] = bisection(IsoboleResponse, @(a, b)C.hand(a, b), aMin, bMin, aMax, bMax);
        end

    end
    
	if strcmp(space, 'log')
        hold on
        loglog(isoboleB, isoboleA);
    else
        hold on
        plot(isoboleB, isoboleA)
    end        
    
    
%% Bliss

    isoboleA = zeros(numel(lambda), 1);
    isoboleB = zeros(numel(lambda), 1);    
    
    
    for i = 1:numel(lambda)

        aMin = lambda(i) * aDoseAxisLimits(1);
        aMax = 2*lambda(i) * aDoseAxisLimits(2);
        bMin =  (1-lambda(i)) * bDoseAxisLimits(1);
        bMax =2*(1-lambda(i)) * bDoseAxisLimits(2);
        
        [isoboleA(i), isoboleB(i)] = bisection(IsoboleResponse, @(a, b)C.bliss(a, b), aMin, bMin, aMax, bMax);


    end

    if strcmp(space, 'log')
        hold on
        loglog(isoboleB, isoboleA);
    else
        hold on
        plot(isoboleB, isoboleA)
	end
    
    
%% HSA

    doseB = C.DrugB.invertDrug(IsoboleResponse);

    doseA = C.DrugA.invertDrug(IsoboleResponse);

    if strcmp(space, 'log')
        hold on
        loglog([doseB; doseB; bDoseAxisLimits(1)], [aDoseAxisLimits(1); doseA; doseA])
        title('Isoboles in log-space')
    else
        hold on
        plot([doseB; doseB; bDoseAxisLimits(1)], [aDoseAxisLimits(1); doseA; doseA])
        title('Isoboles in linear space')
    end


    plot(10^-3, 1.5, 'p')
    legend('Loewe', 'TallaridaLB', 'TallaridaUB', 'Hand', 'Bliss', 'HSA', 'measurement')
        
    %% Add Legend & fomat Plots
    set(findall(gca, 'Type', 'Line'),'LineWidth',2); % Set Line Thickness    
    %legend('Loewe', 'TallaridaLB', 'TallaridaUB', 'Hand', 'Bliss', 'HSA')
    
     
    


end

