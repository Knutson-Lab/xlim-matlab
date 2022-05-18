classdef TcspcLifetimeModel < TcspcModel
    %TCSPCLIFETIMEMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        % Lifetimes - Lifetime components contained in model
        Lifetimes (:,1) LifetimeComponent
        % Anisotropies - Anisotropy components contained in model
        Anisotropies (:,1) LifetimeComponent
        % ColorShift - Value by which IRF should be shifted when convolving to exponential term
        ColorShift (1,1) OptimizationParameter
        % Scatter - Fractional contribution of IRF
        Scatter (1,1) OptimizationParameter
        % BackgroundFraction - Fractional contribution of background curve
        BackgroundFraction (1,1) OptimizationParameter
        % Offset - Number of photons from uniform background
        Offset (1,1) OptimizationParameter
    end

    properties (Access = protected,Hidden)
        % HIDDENLIFETIMES
        HiddenLifetimes
        % HIDDENANISOTROPIES
        HiddenAnisotropies
        % HIDDENCOLORSHIFT
        HiddenColorShift
        % HIDDENSCATTER
        HiddenScatter
        % HIDDENBACKGROUNDFRACTION
        HiddenBackgroundFraction
        % HIDDENOFFSET
        HiddenOffset
    end    

    properties (Dependent, SetAccess = protected)
        ParameterVector (:,1) OptimizationParameter
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

        function model = set.Lifetimes(model,life)
            model.HiddenLifetimes = life;
        end

        function life = get.Lifetimes(model)
            life = model.HiddenLifetimes;
        end

        function model = set.Anisotropies(model,ani)
            model.HiddenAnisotropies = ani;
        end

        function ani = get.Anisotropies(model)
            ani = model.HiddenAnisotropies;
        end

        function model = set.ColorShift(model,col)
            model.HiddenColorShift = col;
        end

        function col = get.ColorShift(model)
            col = model.HiddenColorShift;
        end

        function model = set.Scatter(model,sca)
            model.HiddenScatter = sca;
        end

        function sca = get.Scatter(model)
            sca = model.HiddenScatter;
        end

        function model = set.BackgroundFraction(model,back)
            model.HiddenBackgroundFraction = back;
        end

        function back = get.BackgroundFraction(model)
            back = model.HiddenBackgroundFraction;
        end

        function model = set.Offset(model,off)
            model.HiddenOffset = off;
        end

        function off = get.Offset(model)
            off = model.HiddenOffset;
        end

        function paramv = get.ParameterVector(model)
            param = 
        end    
        
        function trans = simulate(model, irf, bck, ncurves, args)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            arguments
                % MODEL - Model to simulate
                model (1,1) TcspcLifetimeModel
                % IRF - Instrument response function to use
                irf (1,1) TcspcInstrument
                % BCK - Background data to include with simulated model
                bck (1,1) TcspcData
                % NCURVES - Number of curves to simulate from model
                ncurves (1,1) {mustBeNonnegative, mustBeInteger} = 10000
                % args.SimulationMethod- Options for how photons should be generated from curve
                args.SimulationMethod {mustBeMember(args.SimulationMethod,["Photons","Distribution"])} = "Photons"
                % args.ConvolutionMethod- Option for how IRF should be convolved with curves
                args.ConvolutionMethod {mustBeMember(args.ConvolutionMethod,["FFT","Trapezoidal","Simpsons","QuadGK"])} = "FFT"
                % args.InterpolationMethod- Options for how IRF should be interpolated for Gaussian quadrature
                args.InterpolationMethod {mustBeMember(args.InterpolationMethod,["cubicinterp","smoothingspline"])} = "smoothingspline"
            end

            assert(irf.Data == bck);
            values = model.ParameterVector.multiGet("Value");
            c = namedargs2cell(args);
            transdata = sim_tcspc_dks(irf.Data.BinWidth,irf.Data.Data,bck.Data,values,numel(model.Lifetimes),ncurves,c{:});
            trans(ncurves) = TcspcData;
            trans = trans.multiSet("Data",transdata,"BinWidth",irf.Data.BinWidth);
        end
    end
end

