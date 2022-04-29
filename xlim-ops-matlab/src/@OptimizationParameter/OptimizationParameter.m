classdef OptimizationParameter
    % OPTIMIZATIONPARAMETER Abstract container for all different types of optimization routines

    properties
        % VALUE
        Value (1,1)
        % CONSTRAINTS
        Constraints (1,1)
    end
    
    methods
        function param = OptimizationParameter(val,const)
            % OPTIMIZATIONPARAMETER
            arguments
                val (:,1) FitParameter
                const (:,1) BoundaryConstraints
            end
            if(nargin>1)
                noofval = numel(val);
                noofconst = numel(const);
                assert(noofval == noofconst, "Must be equal number of FitParameter objects and BoundaryConstraint objects");

                values = val.multiGet("Value");
                bounds = const.multiGet("ScalarBounds");

                [lb,ub] = bounds.multiGet("LowerBound","UpperBound");
                assert(numel(lb) == numel(ub),"Number of lower bounds and upper bounds must be equal");
                noofbounds = numel(lb);
                
                if(noofbounds~=0)
                    assert(noofbounds == values,"Each bound must correspond to a value")

                    checklower = all(lb<=values);
                    checkupper = all(values<=ub);                

                    assert(checklower && checkupper,"Some Value is not between ScalarBounds");
                end

                param(1:noofval) = OptimizationParameter; 
                for i=1:noofval
                    param(i).Value = val(i);
                    param(i).Constraints = const(i);
                end 
            end
        end
    end

    methods (Abstract)
        multiparam = multiSet(multiparam, args)
        varargout = multiGet(multiparam, prop)
    end
end

