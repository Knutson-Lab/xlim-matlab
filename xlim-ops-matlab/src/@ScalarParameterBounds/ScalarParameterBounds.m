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
                lb (:,1) {mustBeNumeric} = []
                ub (:,1) {mustBeNumeric} = []
            end

            if(nargin>1)
                nooflb = numel(lb);
                noofub = numel(ub);
                noofBounds = max(nooflb,noofub);
                bounds(1:noofBounds) = ScalarParameterBounds;

                if(~isempty(lb) && ~isempty(ub))
                    assert(all(lb <= ub),"Lower bounds must be less than or equal to corresponding upper bounds")
                    if(nooflb>1 && noofub>1)
                        assert(nooflb==noofub,"No of lower bounds and upper bounds provided not equal")   
                    end
                end 

                for i=1:noofBounds
                    if(nooflb>0)
                        if(nooflb>1)
                            bounds(i).LowerBound = lb(i);
                        else
                            bounds(i).LowerBound = lb;
                        end
                    else
                        bounds(i).LowerBound = [];
                    end
                    if(noofub>0)
                        if(noofub>1)
                            bounds(i).UpperBound = ub(i);
                        else
                            bounds(i).UpperBound = ub;
                        end
                    else
                        bounds(i).UpperBound = [];
                    end 
                end
            end
        end

        function varargout = multiGet(multibound,prop)
            arguments
                multibound ScalarParameterBounds
            end
            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["LowerBound","UpperBound"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                varargout{i} = vertcat(multibound.(prop{i}));
            end    
        end
    end
end

