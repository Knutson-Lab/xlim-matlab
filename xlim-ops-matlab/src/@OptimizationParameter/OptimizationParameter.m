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

