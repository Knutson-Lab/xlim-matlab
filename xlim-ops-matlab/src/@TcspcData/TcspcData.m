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

    methods
        function multiTcspc = multiSet(multiTcspc, args)
            % Multiple curve constructor
            arguments
                multiTcspc (:,1) TcspcData
                args.Data (:,:) {mustBeNonnegative,mustBeInteger}
                args.BinWidth (:,1) {mustBePositive,mustBeFinite}
                args.PolarizationAngle (:,1) {mustBeInRange(args.PolarizationAngle,0,90)}
            end

            ntrans = numel(multiTcspc);

            isData = isfield(args,"Data");
            if(isData)
                isData1 = numel(args.Data(1,:)) == 1;
                assert(numel(args.Data(1,:)) == ntrans | isData1, "No of curves and no of Tcspc objects do not match")
            end

            isBinwidth = isfield(args,"BinWidth");
            if(isBinwidth)
                isBinwidth1 = numel(args.BinWidth) == 1;
                assert(numel(args.BinWidth) == ntrans | isBinwidth1, "No of bin widths provided and no of Tcspc objects do not match")
            end

            isPol = isfield(args,"PolarizationAngle");
            if(isPol)
                isPol1 = numel(args.PolarizationAngle) == 1;
                assert(numel(args.PolarizationAngle) == ntrans | isPol1, "No of pol angles provided and no of Tcspc objects do not match")
            end

            for i = 1:ntrans
                if(isData)
                    if (~isData1)
                        multiTcspc(i).HiddenData = args.Data(:,i);
                    else
                        multiTcspc(i).HiddenData = args.Data(:,1);
                    end
                end
                if(isBinwidth)
                    if (~isBinwidth1)
                        multiTcspc(i).HiddenBinWidth = args.BinWidth(i);
                    else
                        multiTcspc(i).HiddenBinWidth = args.BinWidth(1);
                    end
                end
                if(isPol)
                    if (~isPol1)
                         multiTcspc(i).HiddenPolarizationAngle = args.PolarizationAngle(i);
                    else
                         multiTcspc(i).HiddenPolarizationAngle = args.PolarizationAngle(1);
                    end
                end
            end            
        end

        function [varargout] = multiGet(multiTcspc, prop)
            arguments
                multiTcspc (:,1) TcspcData
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Data","BinWidth","PolarizationAngle","NumberOfBins"])}
            end

            assert(numel(prop) <= 4, "Maximum number of properties accepted is 4")

            for i=1:numel(prop)
                switch prop{i}
                    case "Data"  
                        varargout{i} = horzcat(multiTcspc.Data);
                    case "BinWidth"
                        varargout{i} = vertcat(multiTcspc.BinWidth);
                    case "PolarizationAngle"
                        varargout{i} = vertcat(multiTcspc.PolarizationAngle);
                    case "NumberOfBins"
                        varargout{i} = vertcat(multiTcspc.NumberOfBins);
                    otherwise
                        break; %prop validation done in argument block
                end
            end
        end    
    end
end