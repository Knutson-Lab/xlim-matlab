classdef WeightedLifetimeComponent < LifetimeComponent
    % WEIGHTEDLIFETIMECOMPONENT Contains information about the lifetime
    % component and the amplitude value by which it is scaled
    
    properties (Dependent)
        % Amplitude - Amplitude value by which to scale exponential term
        Amplitude (1,1) OptimizationParameter
    end

    properties (Access = protected,Hidden)
        % HIDDENAAMPLITUDE
        HiddenAmplitude
    end
    
    methods
        function param = WeightedLifetimeComponent(args)
            % WEIGHTEDLIFETIMECOMPONENT Construct an instance of
            % WeightedLifetimeComponent with inherited property Lifetime
            % and property Amplitude, which are both FitParameter objects
            arguments
                args.?WeightedLifetimeComponent
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                param.(arc{i}) = arc{i + 1};
            end
        end
        
        function param = set.Amplitude(param,amp)
            param.HiddenAmplitude = amp;
        end

        function amp = get.Amplitude(param)
            amp = param.HiddenAmplitude;
        end

        function multilife = multiSet(multilife, args)
            arguments
                multilife (:,1) LifetimeComponent
                args.Lifetime (:,1) OptimizationParameter
                args.Amplitude (:,1) OptimizationParameter
            end

            nlife = numel(multilife);

            islife = isfield(args,"Lifetime");
            if(islife)
                nlifegiven = size(args.Lifetime);
                isLife1 = nlifegiven == 1;
                assert(nlifegiven == nlife || isLife1, "No of Lifetime parameters and no of lifetimes provided do not match")
            end

            isamp = isfield(args,"Amplitude");
            if(isamp)
                nampgiven = size(args.Amplitude);
                isAmp1 = nampgiven == 1;
                assert(nampgiven == nlife || isAmp1, "No of Amplitude parameters and no of lifetimes provided do not match")
            end

            for i = 1:nlife
                if(islife)
                    if (~isLife1)
                        multilife(i).HiddenLifetime = args.Lifetime(i);
                    else
                        multilife(i).HiddenLifetime = args.Lifetime;
                    end
                end
                if(isamp)
                    if (~isAmp1)
                        multilife(i).HiddenAmplitude = args.Amplitude(i);
                    else
                        multilife(i).HiddenAmplitude = args.Amplitude;
                    end
                end
            end
        end

        function varargout = multiGet(multilife, prop)
            arguments
                multilife (:,1) LifetimeComponent
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Lifetime","Amplitude","AmpTauPair"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                switch prop{i}
                    case "AmpTauPair"
                        [amp,tau] = multiget(multilife,"Amplitude","Lifetime");
                        ntau = numel(tau);
                        varargout(1:2:2 * ntau) = amp;
                        varargout(2:2:2 * ntau) = tau;
                    case "Lifetime"
                        varargout{i} = vertcat(multiTcspc.(prop{i}));
                    otherwise
                        varargout{i} = vertcat(multiTcspc.(prop{i}));
                end
            end
        end 
    end
end

