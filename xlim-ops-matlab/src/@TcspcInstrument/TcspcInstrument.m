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

%         function multiInst = multiSet(multiInst, args)
% 
%             arguments
%                 multiInst (:,1) TcspcInstrument
%                 args.Data (:,:) {mustBeNonnegative,mustBeInteger}
%                 args.ReferenceLifetime
%             end
% 
%         end    
    end 
end


