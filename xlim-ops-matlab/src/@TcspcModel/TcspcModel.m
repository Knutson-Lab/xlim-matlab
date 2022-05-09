classdef (Abstract) TcspcModel
    %TCSPCMODEL Abstract superclass for fitting models to TCSPC curves
    
    properties (Abstract, SetAccess = protected)
        ParameterVector (:,1) OptimizationParameter
        BackgroundFraction (1,1) OptimizationParameter
        Offset (1,1) OptimizationParameter
    end
    
    methods
        function model = TcspcModel(args)
            arguments
                args.?TcspcModel
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                model.(arc{i}) = arc{i + 1};
            end
        end
    end

    methods (Abstract)
        trans = simulate(model);
    end
end

