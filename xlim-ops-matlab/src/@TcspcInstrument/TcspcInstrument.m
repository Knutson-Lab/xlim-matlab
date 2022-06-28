classdef TcspcInstrument
    %TcspcInstrument - Contains data for an IRF acquired using a TCSPC setup
    
    properties (Dependent)
        % DATA - TcspcData acquired on setup
        Data (1,1) TcspcData
        % REFERENCELIFETIME - Lifetime value of known standard acquired on TCSPC system
        ReferenceLifetime (1,1) {mustBePositive, mustBeFinite}
    end

    properties (Access = protected,Hidden)
        % HIDDENDATA
        HiddenData
        % HIDDENREFERENCELIFETIME
        HiddenReferenceLifetime
    end

    methods 
        function irf = TcspcInstrument(args)
            arguments
                args.?TcspcInstrument
            end
            argc = namedargs2cell(args);
            for i = 1:2:numel(argc)
                irf.(argc{i}) = argc{i + 1};
            end
        end

        function irf = set.Data(irf,tcspcdata)
            irf.HiddenData = tcspcdata;
        end
    
        function tcspcdata = get.Data(irf)
            tcspcdata = irf.HiddenData;
        end   

        function irf = set.ReferenceLifetime(irf,life)
            irf.HiddenReferenceLifetime = life;
        end
    
        function life = get.ReferenceLifetime(irf)
            life = irf.HiddenReferenceLifetime;
        end    

        function multiInst = multiSet(multiInst, args)
            % Multiple intrsument response function constructor. Receives
            % arguments of a TcspcInstrument array and name-value pairs for
            % sets of Data (irf data) and ReferenceLifetime to return an
            % array of Tcspc Instrument objects with those properties

            arguments
                multiInst (:,1) TcspcInstrument
                args.Data (:,1) TcspcData
                args.ReferenceLifetime (:,1) {mustBePositive,mustBeFinite}
            end

            nirf = numel(multiInst);

            isData = isfield(args,"Data");
            if(isData)
                ndatacurves = numel(args.Data);
                isData1 = ndatacurves == 1;
                assert(ndatacurves == nirf || isData1,'MATLAB:assertion:failed', "No of irfs and no of TcspcInstrument objects do not match")
            end

            isRefLife = isfield(args,"ReferenceLifetime");
            if(isRefLife)
                nreflife = numel(args.ReferenceLifetime);
                isRefLife1 = nreflife == 1;
                assert(nreflife == nirf || isRefLife1,'MATLAB:assertion:failed', "No of reference lifetimes provided and no of TcspcInstrument objects do not match")
            end

            for i = 1:nirf
                if(isData)
                    if (~isData1)
                        multiInst(i).HiddenData = args.Data(i);
                    else
                        multiInst(i).HiddenData = args.Data;
                    end
                end

                if(isRefLife)
                    if (~isRefLife1)
                        multiInst(i).HiddenReferenceLifetime = args.ReferenceLifetime(i);
                    else
                        multiInst(i).HiddenReferenceLifetime = args.ReferenceLifetime;
                    end
                end
            end             
        end 
        function varargout = multiGet(multiInst, prop)
            %Property parser for single array of TcspcInstrument objects.
            %Receives array f TcspcInstrument and named list of properties
            %to return. Property values parsed and returned as sets of
            %[Data, ReferenceLifetime] in varargout

            arguments
                multiInst (:,1) TcspcInstrument
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Data","ReferenceLifetime"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                varargout{i} = vertcat(multiInst.(prop{i}));
            end
        end         
    end 
end


