classdef TcspcInstrumentModel < TcspcModel
    %TCSPCINSTRUMENTMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        % ScalingFactor - OptimizationParameter that scales Gaussian to match data
        ScalingFactor (1,1) OptimizationParameter
        % Delay - OptimizationParameter that represents time delay between sample excitation and instrument detection
        Delay (1,1) OptimizationParameter
        % FullWidthHalfMax - OptimizationParameter that represents the FWHM of the Gaussian beam
        FullWidthHalfMax (1,1) OptimizationParameter
        % ParameterVector - Bundled OptimizationParameter vector for passage to library routines
        ParameterVector (1,:) OptimizationParameter
    end

    properties (Access = protected,Hidden)
        HiddenScalingFactor
        HiddenDelay
        HiddenFullWidthHalfMax
        HiddenParameterVector
    end  
    
    methods
        function model = TcspcInstrumentModel(args)
            arguments
                args.?TcspcInstrumentModel
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                model.(arc{i}) = arc{i + 1};
            end
        end

        function model = set.ScalingFactor(model,scale)
            model.HiddenScalingFactor = scale;
        end

        function scale = get.ScalingFactor(model)
            scale = model.HiddenScalingFactor;
        end

        function model = set.Delay(model,del)
            model.HiddenDelay = del;
        end

        function del = get.Delay(model)
            del = model.HiddenDelay;
        end

        function model = set.FullWidthHalfMax(model,fwhm)
            model.HiddenFullWidthHalfMax = fwhm;
        end

        function fwhm = get.FullWidthHalfMax(model)
            fwhm = model.HiddenFullWidthHalfMax;
        end

        function model = set.ParameterVector(model,pvec)
            model.HiddenParameterVector = pvec;
        end

        function pvec = get.ParameterVector(model)
            pvec = model.HiddenParameterVector;
        end
        
        function irfs = simulate(model, bck, ncurves, args)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                % MODEL - Instrument Model to simulate
                model (1,1) TcspcInstrumentModel
                % BCK - Background data to include with simulated model
                bck (1,1) TcspcData
                % NCURVES - Number of curves to simulate from model
                ncurves (1,1) {mustBeNonnegative, mustBeInteger} = 1
                % args.SimulationMethod- Options for how photons should be generated from curve
                args.SimulationMethod {mustBeMember(args.SimulationMethod,["Photons","Distribution"])} = "Photons"
                % args.IntegrationOptions- Options for how IRF should be integrated for Gaussian quadrature
                args.IntegrationOptions {mustBeMember(args.IntegrationOptions,["cubicinterp","smoothingspline"])} = "smoothingspline"
            end

            values = model.ParameterVector.multiGet("Value");
            irfdata = sim_tcspc_irfs(bck.Data.BinWidth,bck.Data,values,ncurves,"SimulationMethod",args.SimulationMethod,"IntegrationOptions",args.IntegrationOptions);
            irf_tcspc(ncurves) = TcspcData;
            irf_tcspc = irf_tcspc.multiSet("BinWidth",bck.Data.BinWidth,"Data",irfdata);

            irfs(ncurves) = TcspcInstrument;
            irfs = irfs.multiSet("Data",irf_tcspc);
        end
    end
end

