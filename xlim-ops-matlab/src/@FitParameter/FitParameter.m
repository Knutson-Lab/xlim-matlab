classdef FitParameter
    %FITPARAMETER Structure that holds information about a simulated/fit
    %parameter
    
    properties (Dependent)
        % VALUE - Numeric value of parameter of interest
        Value (1,1) {mustBeNumeric}
        % FIXED - Indicates whether parameter is fixed or free during fit
        Fixed (1,1) logical
        % LINKED - Indicates whether parameter is linked in fitting
        Linked (1,1) logical
    end

    properties (Access = protected,Hidden)
        % HIDDENVALUE
        HiddenValue
        % HIDDENFIXED
        HiddenFixed
        % HIDDENLINKED
        HiddenLinked
    end
    
    methods
        function param = FitParameter(args)
            %FITPARAMETER - Holds the value of a parameter, and whether the
            %parameter is fixed and/or linked
            arguments
                args.?FitParameter
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                param.(arc{i}) = arc{i + 1};
            end
        end
        
        function param = set.Value(param,value)
            param.HiddenValue = value;
        end

        function value = get.Value(param)
            value = param.HiddenValue;
        end

        function param = set.Fixed(param,fix)
            param.HiddenFixed = fix;
        end

        function fix = get.Fixed(param)
            fix = param.HiddenFixed;
        end

        function param = set.Linked(param,link)
            param.HiddenLinked = link;
        end

        function link = get.Linked(param)
            link = param.HiddenLinked;
        end

%         function multiparam = multiSet(multiparam, args)
%             %multiple FitParameter constructor
%             arguments
%                 multiparam (:,1) FitParameter
%                 args.Value (:,1) {mustBeNumeric}
%                 args.Fixed (:,1) logical
%                 args.Linked (:,1) logical
%             end            
% 
%         end

        function varargout = multiGet(multiparam, prop)
            %Property parser for single or array of FitParameter objects
            arguments
                multiparam (:,1) FitParameter
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Value","Fixed","Linked"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                    varargout{i} = vertcat(multiparam.(prop{i}));
            end
        end 
    end
end

