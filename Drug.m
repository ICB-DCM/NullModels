classdef Drug < handle
   %% Summary 
   %
   % handle => only one version and no copies :)
   %
   % The Class supports Hill-Type-Curves:
   %
   %    $f_A(a) = 1 - w_{max} \frac{a^n}{a^n+d^n}$
   %
   %
   %% List of all Properties
   %
   % Name (String)
   % Parameters (Vector) [w_max; d; n];
   %
   % Conc (Vector = Data)
   % Response (Vector = Data)
   % FittingResidual (Double)
   %
   %% List of all methods
   %
   % Drug(newDrugName)                      Constructor
   % addMeasurement(newMeasurement)         
   % fitDrug()                              
   % evaulateObjective(parameters)
   % evaluateDrug(a)
   % invertDrug(e)
   % diffEvaluateDrug(a)
   %
   %
   
properties
    Name % Name of the Drug
    Parameters % Parameters needed for the particular Type
     
    Conc % Vector with the dose-concentrations found in the Data
	Response % Vector imported from Data
    
    FittingResidual % Residual of the fitting;

end
methods
    
    function obj = Drug(newDrugName)
        
	%% Constructor
    %
    % Constructs an empty Drug: 
    %   newDrugName            = Name of Drug (String)
    %   newDrugtype (optional) = Type of Drug (String) , default: 'Hill'
    %
    % Jakob, 23.05.2017
        
        obj.Name = newDrugName;      
        obj.Parameters = [1; 0 ;0];
        
        
        obj.Conc = [];
        obj.Response = [];
        
        obj.FittingResidual = Inf;
    end
    
    %% Needed functionality
     
    
    function this = addMeasurement(this, newConcentration, newResponse)
	%% addMeasurement
    %
    % Adds new Measurements
    %
    % Jakob, 23.05.2017
        
        this.Conc = [this.Conc; newConcentration];
        this.Response = [this.Response; newResponse];
        
    end
    

    
    function this = fitDrug(this)
        %% fitDrugs():
        %
        % Fits the Parameters to the stored Data. Uses lsqnonlin for
        % solving the Lsq-Problem and evaulateObjective for obj. fct. & 
        % gradient.
        % Uses Log-Transformed parameters!
        % Uses 10 Starting points from a Latin Hypercube-Approach as
        % initial values.
        %
        %   The Values are forced to be in the range
        %           0 < w < 1
        %           0 < d_{1/2} < max(Dose)
        %           1 < n < 20
        %
        % Jakob, 24.05.2017-10.07.2017
        
        
        %% maximal and minimal value of the samples:
        
        w_min = 0;
        w_max = 1;
        d_min = min(this.Conc);
        d_max = max(this.Conc); % maximal concentration as maximal half-max
        n_min = 1;
        n_max = 20;
        
        %% Building Initial Values from a Latin Hypercube 
        noSamples = 10; %Number of Samples for multistart optimization;       
        
        initialValues = zeros(3, noSamples);
        w_perm = randperm(noSamples);
        d_perm = randperm(noSamples);
        n_perm = randperm(noSamples);
        
        for i = 1:noSamples
            
            initialValues(1, i) = log((w_max - w_min)*(w_perm(i) - 0.5)/noSamples + w_min); 
            
            initialValues(2, i) = log((d_max - w_min)*(d_perm(i) - 0.5)/noSamples + d_min); 
            
            initialValues(3, i) = log((n_max - n_min)*(n_perm(i) - 0.5)/noSamples + n_min);

        end
        
        
        %% Multistart-Optimization
        
        residuals = zeros (1, noSamples);
        optima = zeros(3, noSamples);
        
        % Set optimizer Options
        %options = optimoptions('lsqnonlin', 'Display', 'off');
        %options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true, 'Display', 'off', 'CheckGradients', true); 
        
        options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true, 'Display', 'off');
                
        %% Call lsqnonlin for all Samples
        for i=1:noSamples
            % test if error was thrown in matlab optimizer
            try 
            [optima(:, i), residuals(i)] = ...
                lsqnonlin(@(Parameters)this.evaulateObjective(exp(Parameters)), initialValues(:, i), ...
                [log(10^-15), log(d_min), 0], [0, log(d_max), log(20)], options);
            catch
                %if error due to abnormailities in the solver occur, just
                %sample new initial values and try again.
                this.fitDrug();
                return
            end
        end
        
        %% Test via AIC/BIC if the Constant 0-Response explains Data better
        
        n = sum(~isnan(this.Response)); % number of Measurements
        
        BIC_Hill = n*log(min(residuals)/n) + 3* log(n); % BIC = n*log(residual/n) + k * log(n)
        
        ZeroResponseResidual = 1-this.Response; % this is the (measurementwise) residual for the constant 1 solution
        
        BIC_ZeroResponse = n*log(norm(ZeroResponseResidual)^2/n); % BIC, k = 0;
        
        
        if BIC_Hill < BIC_ZeroResponse % Hill_Modell ist besser

            [this.FittingResidual, index] = min(residuals);
            this.Parameters = exp(optima(:, index));
            
        else
            this.Parameters = [0; 1; 1];   
            this.FittingResidual = norm(ZeroResponseResidual)^2;
        end

        
    end
    
    function [objectiveFct, jacobian] = evaulateObjective(this, parameters)
        %% evaulateObjective(this, newParameters)
        %
        % evaluates Objective & gradient/Jacobian for the least squares fit
        %
        %   NOTE: computes them for log-Parameter, expects exp(\xi) as
        %   inputs, computes d/(d\xi) f(exp(\xi)) = \napla f(exp(\xi)).*\xi
        %
        % Jakob, 18.6.2017
        
        
        j = 1; % Dimension of the probelem, since 0-Doses-measurements are 
                 % not fitted!
                
        objectiveFct = zeros(length(this.Response), 1); %preallocate memory
        jacobian = zeros(length(this.Response), 3); % is more of a Jacobian
        
        for i = 1:length(this.Response) 
            
            a = this.Conc(i);
            
            if a>0 %tests if dose > 0
                
                w_max = parameters(1);
                d = parameters(2);
                n = parameters(3);
                %x0 = 1;
                
                objectiveFct(j) = 1 - w_max* (a^n / (a^n + d^n)) - this.Response(i);
                                
                jacobian(j, :) = -1* [  w_max * a^n/(a^n + d^n), ...
                      (-w_max*a^n*n*d^n)/((a^n+d^n)^2),...
                      (w_max*n*(a*d)^n*log(a/d))/((a^n+d^n)^2)]; 
                 
                j = j + 1;
            end
        end
        
        % extract the "used" dimensions
        objectiveFct = objectiveFct(1:j-1);
        jacobian = jacobian(1:j-1, :);
        
    end
    
    function effect = evaluateDrug(this, a)
        %% evaluateDrug(a)
        %
        % Evaluates the expected Effect of dose a. This is given by a
        % Hill-Curve:
        %
        % $f_A(a) = 1- w_{max} a^n/(a^n+d^n)$ 
        %
        % where w_max = Parameters(1), d = Parameters(2) und n  = Parameters(3)
        %
        % Jakob, 23.5.2017
        
                        
            w_max = this.Parameters(1);
            d = this.Parameters(2);
            n = this.Parameters(3);
            effect = 1 - w_max* (1 / (1+ (d/a)^n));
            %effect = 1 - w_max* (a^n / (a^n + d^n));
    end
    
    function a = invertDrug(this, e)
    %% invertDrug
    % 
    % Derives the needed Dose for a given effect, depending on the Type of
    % the Drug
    %
    %   $f_A^{-1}(e) = d$
    %
    % Jakob, 23.05.2017

	w_max = this.Parameters(1);
    d = this.Parameters(2);
    n = this.Parameters(3);
    
        x = (1-e);
        
        if w_max == 0 && e==1
            a = 0;
            %disp('Inversion of Zero-Response-Drug with w_max == 0 && e==1 => Drug not invertible')
            return;
        end
    
       % Both tests for higher numerical robustness...  ^^
       if  (w_max - x < 0) || e <= 1-w_max
           a = Inf;
           %disp('Invalid effect in invertDrug: Effect below 1 - w_max');
           return
       elseif e >1 
           a = -Inf;
           %disp('Invalid effect in invertDrug: Effect above 1');
           return
       else
                      
           
           %a = 1/w_max * nthroot(w_max * d^n/(w_max/e - 1), n);
           
           a= d * nthroot(x/(w_max-x), n);
           
       end
    end
    
    function df = diffEvaluateDrug(this, a)
        %% diffEvaluateDrug(a)
        %
        % evaluates the derivative of the Hill Curve with respect to a
        %
        % Jakob, 5.07.2017
        
        if isinf(a)
            df = 0;
            return
        end
        
        w_max = this.Parameters(1);
        d = this.Parameters(2);
        n = this.Parameters(3);
        
        df = -w_max * n*a^(n-1)* d^n/((d^n+a^n)^2);
        
        
        if isnan(df) || isinf(df) 
            % Check for Overflow and use alternative formulation from 
            % f(a) = 1-w \frac{1}{1+(d/a)^n}
            if a==0
                df = 0;
            else
                da = (d/a)^n;
                df = -w_max*n/a *(1/(1+da)^2)*da;
            end
        end
        
    end
    
end
end