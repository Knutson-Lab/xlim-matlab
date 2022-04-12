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
        HiddenPolarizationAngle
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
    end

    methods (Static)
        function arr = multiSet(eps,trans,psi)
            arguments
                eps (:,1) {mustBePositive,mustBeFinite}
                trans (:,:) {mustBeNonnegative,mustBeInteger}
                psi (:,1) {mustBeInRange(psi,0,90)} = 44.7
            end
            assert(nargin > 0);

            % Initialize array
            ntrans = size(trans,2);
            % Check to ensure that dimensions of eps, psi align
            arr(1:ntrans) = TcspcData("BinWidth",eps,"PolarizationAngle",psi);
            for i = 1:ntrans
                arr(i).HiddenData = trans(:,i);
                arr(i).HiddenBinWidth = eps(i);
            end
        end
    end
end