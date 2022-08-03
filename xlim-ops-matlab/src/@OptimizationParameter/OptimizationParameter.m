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
                val (:,1) FitParameter = FitParameter();
                const (:,1) BoundaryConstraints = BoundaryConstraints()
            end
            if(nargin>1)
                noofval = numel(val);
                noofconst = numel(const);
                assert(noofval == noofconst || noofconst == 1,'OptimizationParameter:constructor:ValueAndConstraintsNumber', "Must be equal number of FitParameter objects and BoundaryConstraint objects");

                values = val.multiGet("Value");
                bounds = const.multiGet("ScalarBounds");

                [lb,ub] = bounds.multiGet("LowerBound","UpperBound");
                noofub = numel(ub);
                nooflb = numel(lb);

                if(noofub > 0 && nooflb > 0)
                    if(noofub > 1 && nooflb > 1)
                        % array of lower bounds and upper bounds
                        assert(noofub == nooflb,'OptimizationParameter:constructor:LowerUpperBoundNumber',"Must have same number of lower and upper bounds")
                        assert(noofub == noofval,'OptimizationParameter:constructor:ValuesAndBoundsNumber',"Each bound must correspond to a value")
                    elseif(noofub > 1)
                        % 1 lower bound, array of upper bounds
                        assert(noofub == noofval,'OptimizationParameter:constructor:ValuesAndUpperBoundsNumber',"Each bound must correspond to a value")
                    else
                        % 1 upper bound, array of lower bounds
                        assert(nooflb == noofval,'OptimizationParameter:constructor:ValuesAndLowerBoundsNumber',"Each bound must correspond to a value")
                    end

                    checklower = all(lb<=values);
                    checkupper = all(values<=ub);                

                    assert(checklower && checkupper,'OptimizationParameter:constructor:ValueOutsideBounds',"Some Value is not between ScalarBounds")
                elseif (noofub > 0)
                    % 1 or array of upper bounds, no lower bound
                    assert(noofub == noofval || noofub == 1,'OptimizationParameter:constructor:ValuesAndBoundsNumber',"Each bound must correspond to a value")
                    checkupper = all(values<=ub);                

                    assert(checkupper,"Some Value is not between ScalarBounds")
                elseif (nooflb > 0)
                    % 1 or array of lower bounds, no upper bound
                    assert(nooflb == noofval || nooflb == 1,'OptimizationParameter:constructor:ValuesAndBoundsNumber',"Each bound must correspond to a value")
                    checklower = all(lb<=values);              

                    assert(checklower,'OptimizationParameter:constructor:ValueOutsideBounds',"Some Value is not between ScalarBounds")
                end    

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

