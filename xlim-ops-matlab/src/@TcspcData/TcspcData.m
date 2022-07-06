classdef TcspcData
    % TcspcData - Data container for transient curves collected using time
    % correlated single photon counting. Constructor receives name-value
    % pair arguments for Data, BinWidth, and PolarizationAngle and returns
    % a TcspcData object with those property values

    properties (Dependent)
        % DATA - Holds column vector of photon counts per time channel
        Data (:,1) {mustBeNonnegative,mustBeInteger}
        % BINWIDTH - Time resolution of photon collection bins in nanoseconds
        BinWidth (1,1) {mustBePositive,mustBeFinite}
        % POLARIZATIONANGLE - Angle at which emitted photons from sample were collected
        PolarizationAngle (1,1) {mustBeInRange(PolarizationAngle,0,90)}
    end

    properties (Dependent, SetAccess = protected)
        % NUMBEROFBINS - Number of time bins in the acquired data
        NumberOfBins (1,1) {mustBeNonnegative, mustBeInteger}
    end

    properties (Access = protected,Hidden)
        % HIDDENDATA
        HiddenData
        % HIDDENBINWIDTH
        HiddenBinWidth
        % HIDDENPOLARIZATIONANGLE
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
            % Multiple curve constructor. Receives arguments of a TcspcData
            % array and name-value pairs for sets of Data, BinWidth,
            % Polarizations to return an array of TcspcData objects with
            % those properties
            arguments
                multiTcspc (:,1) TcspcData
                args.Data (:,:) {mustBeNonnegative,mustBeInteger}
                args.BinWidth (:,1) {mustBePositive,mustBeFinite}
                args.PolarizationAngle (:,1) {mustBeInRange(args.PolarizationAngle,0,90)}
            end

            ntrans = numel(multiTcspc);

            isData = isfield(args,"Data");
            if(isData)
                ndatacurves = size(args.Data,2);
                isData1 = ndatacurves == 1;
                assert(ndatacurves == ntrans || isData1,'TcspcData:multiSet:DataSize', "No of curves and no of Tcspc objects do not match")
            end

            isBinwidth = isfield(args,"BinWidth");
            if(isBinwidth)
                nbinwidth = numel(args.BinWidth);
                isBinwidth1 = nbinwidth == 1;
                assert(nbinwidth == ntrans | isBinwidth1,'TcspcData:multiSet:BinWidthSize', "No of bin widths provided and no of Tcspc objects do not match")
            end

            isPol = isfield(args,"PolarizationAngle");
            if(isPol)
                npolang = numel(args.PolarizationAngle);
                isPol1 = npolang == 1;
                assert(npolang == ntrans | isPol1,'TcspcData:multiSet:PolAngSize', "No of pol angles provided and no of Tcspc objects do not match")
            end

            for i = 1:ntrans
                if(isData)
                    if (~isData1)
                        multiTcspc(i).HiddenData = args.Data(:,i);
                    else
                        multiTcspc(i).HiddenData = args.Data;
                    end
                end
                if(isBinwidth)
                    if (~isBinwidth1)
                        multiTcspc(i).HiddenBinWidth = args.BinWidth(i);
                    else
                        multiTcspc(i).HiddenBinWidth = args.BinWidth;
                    end
                end
                if(isPol)
                    if (~isPol1)
                         multiTcspc(i).HiddenPolarizationAngle = args.PolarizationAngle(i);
                    else
                         multiTcspc(i).HiddenPolarizationAngle = args.PolarizationAngle;
                    end
                end
            end            
        end

        function varargout = multiGet(multiTcspc, prop)
            %Property parser for single or array of TcspcData objects.
            %Receives array of TcspcData and named list of properties to
            %return. Property values are parsed and returned as sets of
            %[Data, BinWidth, PolarizationAngle, NumberOfBins] in varargout
            arguments
                multiTcspc (:,1) TcspcData
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Data","BinWidth","PolarizationAngle","NumberOfBins"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                switch prop{i}
                    case "Data" 
                        uniquebins = unique(vertcat(multiTcspc.NumberOfBins));
                        assert(isscalar(uniquebins),'TcspcData:multiGet:BinWidthSizes',"Number of bins not the same in all curves");
                        varargout{i} = horzcat(multiTcspc.Data);
                    otherwise
                        varargout{i} = vertcat(multiTcspc.(prop{i}));
                end
            end
        end    
    end
end