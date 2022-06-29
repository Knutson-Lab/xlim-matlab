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

        function multiparam = multiSet(multiparam, args)
            %multiple FitParameter constructor
            arguments
                multiparam (:,1) FitParameter
                args.Value (:,1) {mustBeNumeric}
                args.Fixed (:,1) logical
                args.Linked (:,1) logical
            end            

            nparam = numel(multiparam);

            isValue = isfield(args,"Value");
            if(isValue)
                nvalue = numel(args.Value);
                isValue1 = nvalue == 1;
                assert(nvalue == nparam | isValue1,'MATLAB:assertion:failed', "No of values provided and no of FitParameter objects do not match")
            end

            isFixed = isfield(args,"Fixed");
            if(isFixed)
                nfixed = numel(args.Fixed);
                isFixed1 = nfixed == 1;
                assert(nfixed == nparam | isFixed1,'MATLAB:assertion:failed', "No of fixed booleans provided and no of FitParameter objects do not match")
            end   

            isLinked = isfield(args,"Linked");
            if(isLinked)
                nlinked = numel(args.Linked);
                isLinked1 = nlinked == 1;
                assert(nlinked == nparam | isLinked1,'MATLAB:assertion:failed', "No of linked booleans provided and no of FitParameter objects do not match")
            end  

            for i = 1:nparam
                if(isValue)
                    if (~isValue1)
                        multiparam(i).HiddenValue = args.Value(i);
                    else
                        multiparam(i).HiddenValue = args.Value;
                    end
                end
                if(isFixed)
                    if (~isFixed1)
                        multiparam(i).HiddenFixed = args.Fixed(i);
                    else
                        multiparam(i).HiddenFixed = args.Fixed;
                    end
                end
                if(isLinked)
                    if (~isLinked1)
                        multiparam(i).HiddenLinked = args.Linked(i);
                    else
                        multiparam(i).HiddenLinked = args.Linked;
                    end
                end
            end
        end

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

