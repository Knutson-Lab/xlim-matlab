classdef TcspcData
    %TCSPCDATA 

    properties (Dependent)
        % DATA - Vector of data
        Data (:,1) {mustBeNonnegative,mustBeInteger}
        % BINWIDTH
        BinWidth (1,1) {mustBePositive,mustBeFinite}
        % POLARIZATIONANGLE
        PolarizationAngle (1,1) {mustBeInRange(PolarizationAngle,0,90)}
    end

    properties (Dependent, SetAccess = protected)
        % Storage for BinWidth
        NumberOfBins (1,1) {mustBePositive,mustBeFinite}
    end

    properties (Access = protected,Hidden)
        % HIDDENDATA
        HiddenData
        % HIDDENBINWIDTH
        HiddenBinWidth
        %HIDDENPOLARIZATIONANGLE
        HiddenPolarizationAngle = 44.7
    end

    methods

        function data = TcspcData(args)
            arguments
                args.?TcspcData
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                data.(arc{i}) = arc{i + 1};
            end
        end

        function data = set.Data(data,trans)
            data.HiddenData = trans;
        end

        function trans = get.Data(data)
            trans = data.HiddenData;
        end

        function data = set.BinWidth(data,binwidth)
            data.HiddenBinWidth = binwidth;
        end

        function binwidth = get.BinWidth(data)
            binwidth = data.HiddenBinWidth;
        end

        function data = set.PolarizationAngle(data,pol)
            data.HiddenPolarizationAngle = pol;
        end

        function polang = get.PolarizationAngle(data)
            polang = data.HiddenPolarizationAngle;
        end

        function bins = get.NumberOfBins(data)
            bins = numel(data.Data);
        end
    end

    methods (Static)

        function multiTcspc = multiSet(args)
            % Multiple curve constructor
            arguments
                args.Data (:,:) {mustBeNonnegative,mustBeInteger}
                args.BinWidth (:,1) {mustBePositive,mustBeFinite}
                args.PolarizationAngle (:,1) {mustBeInRange(args.PolarizationAngle,0,90)}
            end

            ntrans = size(args.Data,2);
            argCell = namedargs2cell(args);
            multiTcspc(1:ntrans) = TcspcData;
            for i = 1:ntrans
                multiTcspc(i).Data = args.Data(:,i);
                multiTcspc(i).BinWidth = args.BinWidth(i);
                multiTcspc(i).PolarizationAngle = args.PolarizationAngle(i);
            end
        end
    end
end