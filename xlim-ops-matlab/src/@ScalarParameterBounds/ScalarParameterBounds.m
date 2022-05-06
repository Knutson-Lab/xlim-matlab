classdef ScalarParameterBounds
    %ScalarParameterBounds - Boundaries around parameter for fitting
    
    properties (SetAccess = protected)
        % LOWERBOUND - Lower bound around parameter
        LowerBound
        % UPPERBOUND - Upper bound around parameter
        UpperBound
    end
    
    methods
%         function bounds = ScalarParameterBounds(lb,ub)
%             % Constructs boundaries for object
%             arguments
%                 lb (:,1) {mustBeNumeric} = -Inf
%                 ub (:,1) {mustBeNumeric} = Inf
%             end
% 
%             if(nargin>1)
%                 noofBounds = numel(lb);
%                 assert(noofBounds == numel(ub),"No of lower bounds and uppder bounds provided not the same");
%                 assert(all(lb <= ub),"Lower bounds must be less than or equal to corresponding upper bounds")
%     
%                 bounds(1:noofBounds) = ScalarParameterBounds;
%                 
%                 for i=1:noofBounds
%                     bounds(i).LowerBound = lb(i);
%                     bounds(i).UpperBound = ub(i);
%                 end
%             end
%         end

function bounds = ScalarParameterBounds(args)
            arguments
                args.LowerBound (:,1) {mustBeNumeric} = []
                args.UpperBound (:,1) {mustBeNumeric} = []
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                bounds.(arc{i}) = arc{i + 1};
            end

            if(nargin>1)
                if(~isempty(args.LowerBound) && ~isempty(args.UpperBound))
                    noofBounds = numel(args.LowerBound);
                    assert(noofBounds == numel(args.UpperBound),"No of lower bounds and uppder bounds provided not the same");
                    assert(all(lb <= ub),"Lower bounds must be less than or equal to corresponding upper bounds")
        
                    bounds(1:noofBounds) = ScalarParameterBounds;
                    
                    for i=1:noofBounds
                        bounds(i).LowerBound = args.LowerBound(i);
                        bounds(i).UpperBound = args.UpperBound(i);
                    end

                elseif(~isempty(args.LowerBound))
                    noofBounds = numel(args.LowerBound);
                    bounds(1:noofBounds) = ScalarParameterBounds;
                    
                    for i=1:noofBounds
                        bounds(i).LowerBound = args.LowerBound(i);
                    end
                else 
                    noofBounds = numel(args.UpperBound);
                    bounds(1:noofBounds) = ScalarParameterBounds;
                    
                    for i=1:noofBounds
                        bounds(i).UpperBound = args.UpperBound(i);
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

