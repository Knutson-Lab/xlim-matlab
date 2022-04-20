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

            if(nargin>1)
                noofBounds = numel(lb);
                assert(noofBounds == numel(ub),"No of lower bounds and uppder bounds provided not the same");
    
                bounds(1:noofBounds) = ScalarParameterBounds;
                
                for i=1:noofBounds
                    assert(lb(i) <= ub(i), '%f must be less than or equal to %f',lb(i),ub(i))
                    bounds(i).LowerBound = lb(i);
                    bounds(i).UpperBound = ub(i);
                end
            end
        end
    end
end

