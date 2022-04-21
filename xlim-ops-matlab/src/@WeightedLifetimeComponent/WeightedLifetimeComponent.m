classdef WeightedLifetimeComponent < LifetimeComponent
    % WEIGHTEDLIFETIMECOMPONENT Contains information about the lifetime
    % component and the amplitude value by which it is scaled
    
    properties (Dependent)
        % Amplitude - Amplitude value by which to scale exponential term
        Amplitude (1,1) FitParameter
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
    end
end

