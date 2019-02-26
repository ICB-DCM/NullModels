function [aDose, bDose] = bisection(responseLevel, responseCurve, aMin, bMin, aMax, bMax)
%% bisection(...)
%
% Computes aDose, bDose, such that response_level = responseCurve(aDose, bDose)
%
% Assumes that aMin/bMin = aMax/bMax, since we want to solve the
% Isobole-bisection along a ray in the dose plane.


    % Check if the initial Values give Dose-Responses on both sides of the
    % responseLevel.

    while responseCurve(aMin, bMin) < responseLevel
        aMin = aMin/2;
        bMin = bMin/2;
        disp('Incorrect initial Values for bisection. Searching for correct values.');
    end

    
    while responseCurve(aMax, bMax) > responseLevel
        aMax = 2*aMax;
        bMax = 2*bMax;
        disp('Incorrect initial Values for bisection. Searching for correct values.');
    end



        aDose = (aMax+aMin)/2;
        bDose = (bMax+bMin)/2;

        while log(aMax)-log(aMin) > 1e-3 || log(bMax)-log(bMin) > 1e-4 % We want to plot in Log-Space => want to be accurate in log space.

            if responseCurve(aDose, bDose) < responseLevel
                aMax = aDose;
                bMax = bDose;
            elseif responseCurve(aDose, bDose) > responseLevel
                aMin = aDose;
                bMin = bDose;
            else
                disp('Error in bisection function')
                keyboard;
            end
            aDose = (aMax+aMin)/2;
            bDose = (bMax+bMin)/2;

        end


end

