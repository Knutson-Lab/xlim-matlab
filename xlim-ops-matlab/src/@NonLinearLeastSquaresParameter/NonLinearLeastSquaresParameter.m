classdef NonLinearLeastSquaresParameter < OptimizationParameter
    % NONLINEARLEASTSQUARESPARAMETER Subclass that is designed for NLLS routines
    methods

%         function param = NonLinearLeastSquaresParameter(val,const)
%             % NONLINEARLEASTSQUARESPARAMETER
%             arguments
%                 val (:,1) FitParameter
%                 const (:,1) BoundaryConstraints = BoundaryConstraints()
%             end
% 
%             param@OptimizationParameter(val,const)
%         end

        function multiparam = multiSet(multiparam, args)

            arguments
                multiparam (:,1) NonLinearLeastSquaresParameter
                args.Value (:,1) FitParameter
                args.Constraint (:,1) BoundaryConstraints
            end

            nparam = numel(multiparam);

            isValue = isfield(args,"Value");
            if(isValue)
                nvalues = length(args.Value);
                isValue1 = nvalues == 1;
                assert(nvalues == nparam || isValue1,'NonLinearLeastSquaresParameter:multiSet:ValueSize', "No of FitParameter values and no of NLLSParameter objects do not match")
            end

            isConstraint = isfield(args,"Constraint");
            if(isConstraint)
                nconsts = length(args.Constraint);
                isconsts1 = nconsts == 1;
                assert(nconsts == nparam || isconsts1,'NonLinearLeastSquaresParameter:multiSet:ConstarintsSize', "No of Boundary Constraints and no of NLLSParameter objects do not match")
            end

            for i = 1:nparam
                if(isValue)
                    if (~isValue1)
                        multiparam(i).Value = args.Value(i);
                    else
                        multiparam(i).Value = args.Value;
                    end
                end
                if(isConstraint)
                    if (~isconsts1)
                        multiparam(i).Constraints = args.Constraint(i);
                    else
                        multiparam(i).Constraints = args.Constraint;
                    end
                end
            end            
        end
        function varargout = multiGet(multiparam, prop)
            arguments
                multiparam (:,1) NonLinearLeastSquaresParameter
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["Value","Constraints"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                    varargout{i} = vertcat(multiparam.(prop{i}));
            end
        end  
    end
end

