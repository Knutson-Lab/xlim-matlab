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
                        for i=1:noofBounds
                            bounds(i).LowerBound = lb(i);
                            bounds(i).UpperBound = ub(i);
                        end
                    elseif(nooflb>1)
                        for i=1:noofBounds
                            bounds(i).LowerBound = lb(i);
                            bounds(i).UpperBound = ub;
                        end
                    else
                        for i=1:noofBounds
                            bounds(i).LowerBound = lb;
                            bounds(i).UpperBound = ub(i);
                        end
                    end
                elseif(~isempty(lb))
                    for i=1:noofBounds
                        bounds(i).LowerBound = lb(i);
                        bounds(i).UpperBound = [];
                    end
                else
                    for i=1:noofBounds
                        bounds(i).LowerBound = [];
                        bounds(i).UpperBound = ub(i);
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

