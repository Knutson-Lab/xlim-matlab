classdef ScalarParameterBounds
    %ScalarParameterBounds - Boundaries around parameter for fitting
    
    properties (SetAccess = protected)
        % LOWERBOUND - Lower bound around parameter
        LowerBound
        % UPPERBOUND - Upper bound around parameter
        UpperBound
    end
    
    methods
        function bounds = ScalarParameterBounds(lb,ub)
            % Constructs boundaries for object
            arguments
                lb (:,1) {mustBeNumeric} = -Inf
                ub (:,1) {mustBeNumeric} = Inf
            end

            noofBounds = numel(lb);


            
            for i=1:noofBounds
                assert(lb(i) <= ub(i), '%f must be less than or equal to %f',lb(i),ub(i))
                bounds(i).LowerBound = lb(i);
                bounds(i).UpperBound = ub(i);
            end
        end
    end
end

