classdef TcspcInstrumentModel < TcspcModel
    %TCSPCINSTRUMENTMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = TcspcInstrumentModel(inputArg1,inputArg2)
            %TCSPCINSTRUMENTMODEL Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function trans = simulate(model, bck, ncurves)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

