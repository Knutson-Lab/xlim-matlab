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
                assert(noofval == noofconst || noofconst == 1, "Must be equal number of FitParameter objects and BoundaryConstraint objects");

                values = val.multiGet("Value");
                bounds = const.multiGet("ScalarBounds");

                [lb,ub] = bounds.multiGet("LowerBound","UpperBound");
                noofbounds = numel(lb);
                
                if(~isempty(lb) && ~isempty(ub))
                    assert(noofbounds == values || noofbounds == 1,"Each bound must correspond to a value")

                    checklower = all(lb<=values);
                    checkupper = all(values<=ub);                

                    assert(checklower && checkupper,"Some Value is not between ScalarBounds");
                elseif(isempty(lb) && ~isempty(ub))
                    assert(noofbounds == values || noofbounds == 1,"Each bound must correspond to a value")

                    checkupper = all(values<=ub);                

                    assert(checkupper,"Some Value is not between ScalarBounds");
                elseif(~isempty(lb) && isempty(ub))
                    assert(noofbounds == values || noofbounds == 1,"Each bound must correspond to a value")

                    checklower = all(lb<=values);              

                    assert(checklower,"Some Value is not between ScalarBounds");
                end

                param(1:noofval) = OptimizationParameter; 
                for i=1:noofval
                    param(i).Value = val(i);
                    if(noofconst==1)
                        param(i).Constraints = const;
                    else 
                        param(i).Constraints = const(i);
                    end    
                end 
            end
        end
    end

    methods (Abstract)
        multiparam = multiSet(multiparam, args)
        varargout = multiGet(multiparam, prop)
    end
end

