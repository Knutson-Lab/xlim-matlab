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

        function multilife = multiSet(multilife, args)
            arguments
                multilife (:,1) LifetimeComponent
                args.Lifetime (:,1) OptimizationParameter
            end

            nlife = numel(multilife);

            islife = isfield(args,"Lifetime");
            if(islife)
                nlifegiven = size(args.Lifetime);
                isLife1 = nlifegiven == 1;
                assert(nlifegiven == nlife || isLife1, "No of Lifetime parameters and no of lifetimes provided do not match")
            end

            for i = 1:nlife
                if(islife)
                    if (~isLife1)
                        multilife(i).HiddenLifetime = args.Lifetime(i);
                    else
                        multilife(i).HiddenLifetime = args.Lifetime;
                    end
                end
            end
        end 

        function varargout = multiGet(multilife, prop)
            arguments
                multilife (:,1) LifetimeComponent
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Lifetime","AmpTauPair"])}
            end

            assert(numel(prop)==1,"LifetimeComponent class has no amplitude")

            switch prop
                case "AmpTauPair" 
                    varargout = multiGet(multilife, "Lifetime");
                otherwise
                    varargout = vertcat(multiTcspc.(prop));
            end
        end 
    end
end

