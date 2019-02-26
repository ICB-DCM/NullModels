classdef Data < handle
    %% Summary
    %
    %   Contains all the Data (= a array of CellLines) from a given data
    %   set. Here single-Drug and combination data are given seperately.
    %
    %% List of all properties:
    %
    %   CellLines:          Array of all cellLines
    %
    %% List of all methods:
    %
    % Data(singleDrugData, combinationData)             Constructor
    %
    % readInDrugs(...)
    % readInCombinations(...)
    % getCellLinePosition(...)
    %
    
    
    properties
        CellLines % cell array containing all Cell Lines in the Data
    end
    
    methods
        
        function obj = Data(singleDrugData, combinationData)
        %% Data(...)
        %
        % Constructor, takes the directory to the files of the single Agent
        % Data and the combination Data. 
        %

        obj.CellLines = cell(0); % initialize the CellLines object with zeros;
        
        if nargin==0 % Create empty Data
            return
        end
        
        %% Read in the single agent data

        obj.readInDrugs(singleDrugData);

        
        %% read in combination data
         
        obj.readInCombinations(combinationData);
            
        end
        
        function this = readInDrugs(this, singleDrugData)
            %% readInDrugs(singleDrugData)
            %
            %   Reads in the single drug data from the file, where
            %   singleDrugData gives the directory of that file. Uses
            %   the viability as readout and fits the hill coeffitients.
            %   The Data Format is in accordance to O'Neil et al. 2016
            %
            %
            % Jakob, 17.4.2018
            
            [num, txt, ~] = xlsread(singleDrugData); % num contains numeric infos, txt the needed string informations
                                                                                % IMPORTANT: txt has one line more for the column-names 
                                                                                %       => indexes are different between num and txt
            
            [noDataPoints, ~]  = size(num);
            
            i = 1;
            
            while i < noDataPoints
                % While you have data from the same Cell Line for the same
                % drug
                
                 j = 1;
                 
                 cellLineName = string(txt{i+1, 2}); % Name of the current cell line
                 drugName = string(txt{i+1, 3}); % Name of the current drug
                 
                 while  i+j <= noDataPoints && strcmp(cellLineName, string(txt{i+j + 1, 2}))  && strcmp(drugName, string(txt{i+j+1, 3}) )
                     j = j+1;
                 end
                 
                 % if we take viability as the measurement of the drug
                 % response
                 
                 drugDose = repmat(num(i:i+j-1, 4), 6, 1); % repmat for the 6 replicates of each experiment
                 drugResponse = reshape(num(i:i+j-1, 5:10), [], 1); % reshape to get the responses as vector
                 
                 %filter out the missing replicates
                 drugDose = drugDose(~isnan(drugResponse));
                 drugResponse = drugResponse(~isnan(drugResponse));
                 
                 % we take mu/muMax as readout
                 %drugDose = num(i:i+j-1, 4); %
                 %drugResponse = num(i:i+j-1, 11);
                 
                 pos = this.getCellLinePosition(cellLineName);
                 
                 if pos==0 % if the cell line is new
                     
                     pos = length(this.CellLines)  + 1 ;
                     this.CellLines{pos} = CellLine(cellLineName);
                     
                 end
                 
                 % Add Measurements
                 this.CellLines{pos}.addDrugMeasurement(drugName, drugDose, drugResponse);
                 
                 i = i+j;
            end
            
            % fit Drugs
            for i = 1:length(this.CellLines)
                this.CellLines{i}.fitDrugs();
            end
            
        end
                
        
        function this = readInCombinations(this, combinationData)
            %% readInCombinations(...)
            %
            %   Reads in the Combinations from the data file (directory) given in
            %   combinationData. Uses mu/muMax as readout.
            %
            % Jakob 19.4.2018
            
            [num, txt, ~] = xlsread(combinationData); % num contains numeric infos, txt the needed string informations
                                                                                % IMPORTANT: txt has one line more for the column-names 
                                                                                %       => indexes are different between num and txt
            
            [noDataPoints, ~]  = size(num);
            
            i = 1;
            
            while i < noDataPoints
            
                j = 1;
                % i+j = index in the numeric- part
                 
                
                
                cellLineName = string(txt{i+1, 2}); % Name of the current cell line
                drugAName = string(txt{i+1, 3}); % Name of the drug A
                drugBName = string(txt{i+1, 5}); % Name of the drug B
                
                
                % we extract which block does belong to the same
                % combination (in the same cell line...)
                
                 while  i+j <= noDataPoints && strcmp(cellLineName, string(txt{i+j + 1, 2}))  && ...
                         strcmp(drugAName, string(txt{i+j+1, 3}) ) && strcmp(drugBName, string(txt{i+j+1, 5}))
                     j = j+1;
                 end
                 
                 
                 drugADose = num(i:i+j-1, 4); % Dose of Drug A
                 drugBDose = num(i:i+j-1, 6); % Dose of Drug B
                 
                 %drugResponse = reshape(num(i:i+j-1, 8:11), [], 1); % Read out.
                 drugResponse = num(i:i+j-1, 8:11); % Read out
                 
                 
                 
                % search for cell Line and add combination
                
                pos = this.getCellLinePosition(cellLineName);
                 
                 if pos==0 % if the cell line is new
                     
                     disp('Error: Cell Line in combination file unknown');
                     disp(cellLineName);
                     
                 end
                 
                 %add Combination Measurement
                this.CellLines{pos}.addCombinationMeasurement(drugAName, drugBName, drugADose, drugBDose, drugResponse);
                
                i = i+j;
            end
            
            
        end
        
                
        function pos = getCellLinePosition(this, testName)
            %% getCellLinePosition(testName)
            %
            % Searches for the Position of the CellLines with Name testName
            % (String) in the CellLines-Cell-Array. Returns Index or 0, if
            % CellLine is not found.
            %
            %
            % Jakob, 17.04.2018
            
            pos = 0;
            for i = 1: length(this.CellLines)
                if strcmp(this.CellLines{i}.Name, testName)
                    pos = i;
                    return;
                end
            end        
        end
        
    end
    
end

