classdef LifetimeComponent
    % LIFETIMECOMPONENT Contains information about lifetime component
    
    properties (Dependent)
        % LIFETIME - Parameter in the denominator of the lifetime term
        Lifetime (1,1) OptimizationParameter
    end

    properties (Access = protected,Hidden)
        % HIDDENLIFETIME
        HiddenLifetime
    end
    
    methods
        function param = LifetimeComponent(args)
            % LIFETIMECOMPONENT Construct an instance of LifetimeComponent
            %with the property Lifetime which is a FitParameter
            arguments
                args.?LifetimeComponent
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                param.(arc{i}) = arc{i + 1};
            end
        end
        
        function param = set.Lifetime(param,life)
            param.HiddenLifetime = life;
        end

        function life = get.Lifetime(param)
            life = param.HiddenLifetime;
        end
    end
end

