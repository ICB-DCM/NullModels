classdef CellLine < handle
        %% Summary 
        %
        %   CellLine contains all the data (= Name, the drugs with the
        %   parameters fitted for this particular cell line and the
        %   combinations with the synergy-coeffitients for this oarticular
        %   cell line) for one cell line
        %
        %% List of all properties:
        %
        % Name 
        % Drugs 
        % Combinations 
        %
        %% List of all Methods:
        %        
        % CellLine(name)                        Constructor
        % addDrugMeasurement(...)
        % fitDrugs() 
        % addCombinationMeasurement()
        %
        %       Auxiliary functions:
        % 
        % getDrugPosition(testName)
        % getCombinationPosition(...)
        %
        
    properties
    
        Name % Name of the cell line
        
        Drugs % Cell-Array of Drugs used in the Experiments
        Combinations % Array/Cell-Array of Combinations found in the Data
        
    end
    
    methods       
        
        
        function obj = CellLine(name)
            %% CellLine(name)
            %
            % Creates a new cell line object
            % Jakob, 17.4.2018
            
            obj.Name = name;
            obj.Drugs = cell(0);
            obj.Combinations = cell(0);           
            
        end
        
        function this = addDrugMeasurement(this, drugName, drugDose, drugResponse)
            %% addDrugMeasurement(...)
            %
            % adds a measurement to the corresponding drug. If there is no
            % data for this Drug in the Cell Line, then the drug is created
            %
            % Jakob, 17.4.2018
            
            pos = this.getDrugPosition(drugName);
            if pos == 0
                
                pos = length(this.Drugs)+1;
                this.Drugs{pos} = Drug(drugName);
                
            end
            
            this.Drugs{pos}.addMeasurement(drugDose, drugResponse);
            
        end
        
        function this = fitDrugs(this)
            %% Fits all Drugs
            %  
            %   fit all Drugs in one CellLine
            % Jakob, 17.4.2018
            
            for i = 1:length(this.Drugs)                
                this.Drugs{i}.fitDrug();
            end    
            
        end
        
        function this = addCombinationMeasurement(this, drugAName, drugBName, ConcA, ConcB, Response)
            %% addCombinationsMeasurement(...)
            %
            %   Adds measurement(s) to the particular combination 
            %
            
            
            pos = this.getCombinationPosition(drugAName, drugBName);
            
            if 0 == pos % test if combination is new, if yes, create this combination
                
                drugA = this.Drugs{this.getDrugPosition(drugAName)};
                drugB = this.Drugs{this.getDrugPosition(drugBName)};
                pos = length(this.Combinations) +1;
                
                this.Combinations{pos} = Combination(drugA, drugB);
            end
            
            % Check in which order the drugs are stored and add Measurement
            % correspondingly
            
            if 1 == this.Combinations{pos}.equalDrugs(drugAName, drugBName)
                    this.Combinations{pos}.addMeasurement(ConcA, ConcB, Response); % add Measurements
            else
                this.Combinations{pos}.addMeasurement(ConcB, ConcA, Response);
            end
            
        end
        
        %% auxiliary functions 
        
        function pos = getDrugPosition(this, testName)
            %% getDrugPosition(testName)
            %
            % Searches for the Position of the Drug with Name testName
            % (String) in the Drugs-Cellarray. Returns Index or 0, if Drug
            % is not found.
            %
            %   getDrugPosition(testName) = i <=> Drugs{i}.Name = testName (or 0!!)
            %
            %
            % Jakob, 23.05.2017
            
            pos = 0;
            for i = 1: length(this.Drugs)
                if strcmp(this.Drugs{i}.Name, testName)
                    pos = i;
                    return;
                end
            end        
        end
        
        
        function pos = getCombinationPosition(this, drugAName, drugBName)
            %% getCombinationPosition(drugAName, drugBName)
            %
            % Function Checks, if there is a Combination with those two 
            % Drugs involved and gives the position. Returns 0 otherwise.
            %
            % Jakob, 17.4.2018
            
            pos = 0;
            for i = 1:length(this.Combinations)
                % Test if right Drugs are involved
                if this.Combinations{i}.equalDrugs(drugAName, drugBName)                    
                    pos = i;
                    return;
                end
            end
        end
        
        %% This functions evaluate the corresponding synergy for their respectiv null modell:
        
        function this = evaluateSynergyLoewe(this)
            for i = 1:length(this.Combinations)
                this.Combinations{i}.evaluateSynergyLoewe();
            end            
        end
        
        function this = evaluateSynergyBliss(this)
            for i = 1:length(this.Combinations)
                this.Combinations{i}.evaluateSynergyBliss();
            end
        end
        
        function this = evaluateSynergyHand(this)
            for i = 1:length(this.Combinations)
                this.Combinations{i}.evaluateSynergyHand();
            end
        end
        
        function this = evaluateSynergyHSA(this)
            for i = 1:length(this.Combinations)
                this.Combinations{i}.evaluateSynergyHSA();
            end
        end
        
        function this = evaluateSynergyTallarida(this)
            for i = 1:length(this.Combinations)
                this.Combinations{i}.evaluateSynergyTallarida();
            end
        end
        
        
    end
    
end