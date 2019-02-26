classdef Combination < handle
    %% Summary
    %
    % An Object of type Combination stores all samples corresponding to a
    % certain Combination.
    %
    %
    %%       List of all Properties
    %
    %        DrugA (Drug)
    %        DrugB (Drug)
    %
    %        ConcA  (Vector)
    %        ConcB  (Vector)
    %        NumberOfCells (Vector)
    %        Response (Matrix)
    %
    %        LoeweIndex (Vector)
    %        LoeweRelativeIndex (Vector)
    %        BlissIndex (Vector)
    %        HandIndex (Vector)
    %        HSAIndex  (Vector)
    %        TallaridaIndex (2xNOMeasurement Vector... [LB, UB])
    %
    %       LoewePrediction
    %       BlissPrediction
    %       HandPrediction
    %       HSAPrediction
    %       TallaridaPrediction (LB, UB)
    %
    %
    %%       List of all methods
    %
    %        Combination(newDrugA, newDrugB)         %Constructor
    %        equalDrugs(drugAName, drugBName)
    %        addMeasurement(this, newMeasurement)
    
    %        evaluateSynergy()
    %        evaluateSynergyLoewe()
    %        evaluateSynergyBliss()
    %        evaluateSynergyHand()
    %        evaluateSynergyHSA()
    %        evaluateSynergyTallarida()
    %
    %
    %        hand()
    %        classicalLoewe()
    %        bliss()
    %        tallaridaAB()
    %        tallaridaBA()
    %
    %        equalDrugs()
    
    
    properties
        DrugA %both Drugs 
        DrugB %
        
        
        %% The following block is imported from the excel-files
        ConcA
        ConcB
        Response
        
        %% The Combination Indizes of all used null models
        
        LoeweIndex % Vector containing the Loewe indizes of all Combinations
        LoeweRelativeIndex % Vector relative Loewe Indices (= prediction / measurement)
        
        BlissIndex % Vector containing the Bliss indizes of all Combinations
        HandIndex % Vector containing the Hand indizes of all Combinations
        HSAIndex  % Vector containing the HSA indizes of all Combinations
        TallaridaIndex % LB & UB Containing the lower and the upper possible value of the Tallarida CI
        
        %% The predicitons of all used null models
        
		LoewePrediction
        BlissPrediction
    	HandPrediction
    	HSAPrediction
        TallaridaPrediction % LB & UB Containing the lower and the upper possible value of the Tallarida Prediction

    end
    
    methods
        
        
        function obj = Combination(newDrugA, newDrugB)
            
        %% Constructor:
        %
        % Constructs new empty Combination
        %
        %   newDrugA /newDrugB = DrugA/DrugB:   Object of class Drug
        %   newSynergyType = SynergyType of the Combination (String, 'Loewe'/'Bliss'); 
        %
        % Jakob, 23.05.2017
            
            obj.DrugA = newDrugA;
            obj.DrugB = newDrugB;
            
            % Empty, has to be imported from the data
            
            obj.ConcA = [];
            obj.ConcB = [];
            obj.Response = [];
            
            obj.LoeweIndex = [];
            obj.LoeweRelativeIndex = [];
            obj.BlissIndex = [];
            obj.HandIndex = [];
            obj.HSAIndex = [];
            obj.TallaridaIndex = [];
            
            obj.LoewePrediction = [];
            obj.BlissPrediction = [];
            obj.HandPrediction = [];
            obj.HSAPrediction = [];
            obj.TallaridaPrediction = [];
            
        end
        
        %% Needed Functionality:
               
        function this = addMeasurement(this, newConcA, newConcB, newResponse) 
        %% addMeasurment
        %
        % Adds new Measurements from the Vector newMeasurement, where
        %
        %   newMeasurement(1) = ConcA
        %   newMeasurement(2) = ConcB
        %   newMeasurement(3) = Response (=Mean_Ch2)
        %
        %   The Loewe and Bliss index remains [] !!!
        %
        %   Jakob, 23.5.2017
            
            this.ConcA = [this.ConcA; newConcA];
            this.ConcB = [this.ConcB; newConcB];
            this.Response = [this.Response; newResponse];
            
        end    
        
        %% Synergy
        %
        %   The Synergy-Functions evaluate the Synergy-Indizes and the
        %   predictions for the corresponding null-model for all
        %   measurement. For a given Dose-Pair the (if multiple
        %   measurements are available... average) Combination Index and
        %   the Response for the Dose Pairs is computed and updated in the
        %   corresponding Response and CI Vector.
        %
        %
        
        
        function this = evaluateSynergy(this)
        %% evaluateSynergs
        %
        % Computes the Synergies (= LoeweIndex / = BlissIndex etc. ) of the given
        % Combination. Takes no input values
        %
        % 19.4.2018 Jakob

            this.evaluateSynergyLoewe();
            this.evaluateSynergyBliss();
            this.evaluateSynergyHand();  
            this.evaluateSynergyHSA();
            this.evaluateSynergyTallarida();
            
        end
        
        
        function this = evaluateSynergyLoewe(this)
        %% evaluateLoewe
        %
        % Computes the Synergies = LoeweIndex of the given Combination
        %
        % 19.4.2018 Jakob

            [noConcs, noReplicates] = size(this.Response); % number of Concentrations/Replicates
        
            this.LoeweIndex = zeros(noConcs, noReplicates);
            this.LoeweRelativeIndex = zeros(noConcs, noReplicates);
                        
            for i = 1:noConcs
                
                %% Loewe
                for j = 1:noReplicates
                    
                    if isnan(this.Response(i, j)) % if there was no corresponding replicate for this Combination
                        this.LoeweIndex(i, j) = nan;
                    else
                        
                        this.LoeweIndex(i, j) = ...
                        (this.ConcA(i))/this.DrugA.invertDrug(this.Response(i, j)) + ...
                        (this.ConcB(i))/this.DrugB.invertDrug(this.Response(i, j));                    
                    
                    end
                end
                
                % Store the Loewe Prediction
                this.LoewePrediction(i) = this.classicalLoewe(this.ConcA(i), this.ConcB(i));
                
                % Also compute the relative Loewe Index
                this.LoeweRelativeIndex(i, :) = (1- this.LoewePrediction(i)) ./ (1 - this.Response(i, :));
                
            end
            
            % compute the average CI for a combination
            this.LoeweIndex = nansum(this.LoeweIndex, 2) ./ sum(~isnan(this.LoeweIndex), 2);
            this.LoeweRelativeIndex = nansum(this.LoeweRelativeIndex, 2) ./ sum(~isnan(this.LoeweRelativeIndex), 2);
            
        end
        
        function this = evaluateSynergyBliss(this)
        %% evaluateSynergsBliss
        %
        % Computes the  BlissIndex of the given combination.
        %
        % 19.4.2018 Jakob

            [noConcs, noReplicates] = size(this.Response); % number of Concentrations/Replicates
                    
            this.BlissIndex = zeros(noConcs, noReplicates);
            
            for i = 1:noConcs
 
                %% Bliss 
                
                 effectOfa = this.DrugA.evaluateDrug(this.ConcA(i));
                 effectOfb = this.DrugB.evaluateDrug(this.ConcB(i));
                
                 % eA & eB have 1 as max and 0 as min. To evaluate the
                 % "probability of binding" in Bliss, we use the Bliss
                 % formula for (1-eA) and (1-eB)a and transform back to the
                 % "1 is Max, 0 is Min"-World:
                 %
                 %       prediction = 1-[(1-eA)+(1-eB) - (1-eA)*(1-eB)]
                 
                 this.BlissPrediction(i) = effectOfa*effectOfb;
                 
                 this.BlissIndex(i, :) = (1- effectOfa*effectOfb) ./ (1 - this.Response(i, :));
                 
            end
            

            
            % compute the average CI for a combination
            this.BlissIndex = nansum(this.BlissIndex, 2) ./ sum(~isnan(this.BlissIndex), 2);
            
        end
        
        function this = evaluateSynergyHand(this)
        %% evaluateSynergsHand
        %
        % Computes the  HandIndex of the given combination.
        %
        % 19.4.2018 Jakob
        
            [noConcs, noReplicates] = size(this.Response); % number of Concentrations/Replicates

            this.HandIndex = zeros(noConcs, noReplicates);
            
            for i = 1:noConcs
                
                this.HandPrediction(i) = this.hand(this.ConcA(i), this.ConcB(i)) ;
                this.HandIndex(i, :) = (1- this.HandPrediction(i)) ./(1 - this.Response(i, :));
            
            end

            % Compute Average CI for all Replicates
            this.HandIndex = nansum(this.HandIndex, 2) ./ sum(~isnan(this.HandIndex), 2);
            
        end
        
        function this = evaluateSynergyHSA(this)
        %% evaluateSynergyHSA()
        %
        %   Evaluates the Synergy for the Highest-Single-Agent reference
        %   Model
        %
        %
        %   Jakob, 24.4.2018
            
            [noConcs, noReplicates] = size(this.Response); % number of Concentrations/Replicates

            this.HSAIndex = zeros(noConcs, noReplicates);
            
            for i = 1:noConcs % Loop over al Dose Pairs (Vecorized => no loop over replicates needed.)
                
                this.HSAPrediction(i) = min( this.DrugA.evaluateDrug(this.ConcA(i)), ...
                                                                this.DrugB.evaluateDrug(this.ConcB(i)) );
                
                this.HSAIndex(i, :) = (1- this.HandPrediction(i)) ./ (1 - this.Response(i, :));
            
            end
            
            % Compute Average CI for all Replicates
            this.HSAIndex = nansum(this.HSAIndex, 2) ./ sum(~isnan(this.HSAIndex), 2);
            
        end
        
        function this =evaluateSynergyTallarida(this)
        %% evaluateSynergyTallarida()
        %
        % evaluates the (average) Combination Indizes (LB & UB) and the
        % predictions (LB & UB) for the Tallarida Null-Model.
        %
        %
                
            [noConcs, noReplicates] = size(this.Response); % number of Concentrations/Replicates

            IndexLB = zeros(noConcs, noReplicates);
            IndexUB = zeros(noConcs, noReplicates);
            
            for i = 1:noConcs % Loop over al Dose Pairs (Vecorized => no loop over replicates needed.)
                
                predictionAB = this.tallaridaAB(this.ConcA(i), this.ConcB(i));
                
                predictionBA = this.tallaridaBA(this.ConcA(i), this.ConcB(i));
                
                this.TallaridaPrediction(i, 1) = max(predictionAB, predictionBA); % the lower Bound
                this.TallaridaPrediction(i, 2) = min(predictionAB, predictionBA); % the upper Bound 
                
                % (since higher values correspond to lower Drug-Responses the min and the max is unintuitively placed!)
                
                IndexLB(i, :) = (1- this.TallaridaPrediction(i, 1)) ./ (1 - this.Response(i, :));
                IndexUB(i, :) = (1- this.TallaridaPrediction(i, 2)) ./ (1 - this.Response(i, :));
                % The Lower and the upper bound for each Response.
                
            end
            
            % Compute Average CI for all Replicates
            this.TallaridaIndex(:, 1) = nansum(IndexLB, 2) ./ sum(~isnan(IndexLB), 2);
            this.TallaridaIndex(:, 2) = nansum(IndexUB, 2) ./ sum(~isnan(IndexUB), 2);
            
        end
        
        %% Dose Response Curves
        
        
        function response = hand(this, aDose, bDose, x0)
        %% hand(aDose, bDose)
        %
        % Computes the hand curve for both drugs and the dosepair
        % (a, b) with an initial value of x0.
        %
        % Jakob 8.9.2017
        
        if aDose + bDose == 0 % Zero Dose => Zero Response
            response = 1;
            return
        end
        
        %% Coordinate Trafo:
        %
        % Since the maximal effect is 0 and the minimal 1 for the readout
        % in O'Neil et al. we coose a coordinate transformation
        %
        %               y = 1-x <=> x = 1-y
        %
        % and compute in the y-space. Thereby we can exploit the
        % theoretical property x<1 => y>0 and set this as an inequality for
        % the numerical ODE solver. This massively improved the stability
        % of our computaitons.
       
        
        if nargin==3
            %x0 = 1-1e-4;
            y0 = 1e-7;
        else
            y0 = 1-x0;
        end
        
        %
        % y = 1-x    <=> x = 1-y
        % \dot{y} = - \dot{x}  
        
        %ODE in y-space
        f = @(t, y) - (aDose/(aDose + bDose) *this.DrugA.diffEvaluateDrug(this.DrugA.invertDrug(1-y)) ... 
              +bDose/(aDose + bDose) *this.DrugB.diffEvaluateDrug(this.DrugB.invertDrug(1-y)));
    
          % Set ODE Solver tolerances and Settings
            options = odeset('AbsTol',1e-8, 'RelTol', 1e-5, 'NonNegative', 1);
            
            %[t, x] = ode15s(f, [0, aDose+bDose], x0, options);
            [t, y] = ode15s(f, [0, aDose+bDose], y0, options);
            
            %response = x(end);
            response = 1-y(end);
                       
        end
         
        function response = classicalLoewe(this, aDose, bDose)
        %% classicalLoewe(aDose, bDose)
        %
        % This function gives the classical Loewe solution to the doses
        % aDose and bDose via solving the equation
        %
        %  $a/f_A^{-1}(x) + b/f_B^{-1}(x) = 1$
        %
        % for x using bisektion.
        %
        % Jakob 17.10.17
        
           % Check if the Effect can be "translated", if not, take the
           % Effect of the drug, that is not "limmiting"
        
           if this.DrugA.evaluateDrug(aDose) < 1 - this.DrugB.Parameters(1)
                response = this.DrugA.evaluateDrug(aDose);
                return
           elseif this.DrugB.evaluateDrug(bDose) < 1 - this.DrugA.Parameters(1)
               response = this.DrugB.evaluateDrug(bDose);
               return
           end
           
           %%
           
            
            response_ub = 1; % upper bound on the Response, here the right hand side of the eqn. is >1
            response_lb = 1 - min(this.DrugA.Parameters(1), this.DrugB.Parameters(1)); % lower bound on the Response, here the right hand side is <1
            
            response = 1/2 * (response_lb + response_ub); % To avoid numerical problems of the condition in the while loop is never fulfilled, which also means that 
                                                                                                %one of the drugs basical has no effect at all.
            

            % Solve the Loewe Eqn. via bisection
            
            while (abs(response_ub - response_lb) > 1e-8) && abs(aDose / this.DrugA.invertDrug(response) + bDose /this.DrugB.invertDrug(response) - 1)>0.001 
                response = 1/2 * (response_lb + response_ub);
                
                if aDose/this.DrugA.invertDrug(response) + bDose/this.DrugB.invertDrug(response) > 1 % caution, the right hand site is monotone increasing w.r.to x !!!
                    response_ub = response;
                else
                    response_lb = response;
                end
                
            end
            

            
        end
        
        function response = bliss(this, aDose, bDose)
        %% Bliss(aDose, bDose)
        %
        %   Evaluates the Bliss-Prediction for a given Dose-Pair
        %
        % Jakob, 8.5.2018
        
        response = this.DrugA.evaluateDrug(aDose) * this.DrugB.evaluateDrug(bDose);
        
        
        end
        
        function response = tallaridaAB(this, aDose, bDose)
        %% tallaridaAB(aDose, bDose)
        %
        % This gives the combined effect as defined by tallarida:
        %
        %   $f_{AB}(a, b) := f_A( a + f_A^{-1}(f_B(b)))$
        %
        % for defining the effect the other way round, see the function
        % tallaridaBA
        %
        % Jakob, 30.04.2018
            
            equipotentDoseOfb = this.DrugA.invertDrug(this.DrugB.evaluateDrug(bDose)); % the equipotent Dose of b
            
            if isinf(equipotentDoseOfb)
                
                if this.DrugB.evaluateDrug(bDose) < 1-this.DrugA.Parameters(1) % f_B(b) exceeds the domain of f_A
                    response = this.DrugB.evaluateDrug(bDose);
                    return
                end
                disp('Error in function tallaridaAB()')
                response = nan;
            else
                response = this.DrugA.evaluateDrug(aDose + equipotentDoseOfb);
            end
        end
        
        function response = tallaridaBA(this, aDose, bDose)
            %% tallaridaBA(aDose, bDose)
            %
            % this is equivalent to tallaridaAB(aDose, bDose), just with the
            % role of both drugs exchanged.
            %
            %   $f_{AB}(a, b) := f_B( b + f_B^{-1}(f_A(a)))$
            %
            % Jakob, 30.04.2018
            
            equipotentDoseOfa = this.DrugB.invertDrug(this.DrugA.evaluateDrug(aDose));
            
            if isinf(equipotentDoseOfa) 
                
                if this.DrugA.evaluateDrug(aDose) < 1-this.DrugB.Parameters(1) % f_A(a) exceeds the domain of f_B
                    response = this.DrugA.evaluateDrug(aDose); 
                    return;
                end
                    disp('Error in function tallaridaBA()')
                    response = nan;
            else
                response = this.DrugB.evaluateDrug(bDose + equipotentDoseOfa);
            end
        
        end
        
        %% auxiliary functions
        
        function res = equalDrugs(this, drugAName, drugBName)
        %% equalDrugs(...)
        %
        % Not the right Drugs =>         =  0 
        % right Drugs, right Order =>  =  1
        % right Drugs wrong Order=> =-1
        %
        % Jakob, 19.4.2018
        
        
            if strcmp(this.DrugA.Name, drugAName) && strcmp(this.DrugB.Name, drugBName)
                res = 1;
                return
            elseif strcmp(this.DrugA.Name, drugBName) && strcmp(this.DrugB.Name, drugAName)
                res = -1;
                return
            end
                res = 0; 
        end
 
    end
    
end

