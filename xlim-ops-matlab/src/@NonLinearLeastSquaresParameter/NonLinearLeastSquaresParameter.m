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

%         function param = NonLinearLeastSquaresParameter(val,const)
%             % NONLINEARLEASTSQUARESPARAMETER
%             arguments
%                 val (:,1) FitParameter
%                 const (:,1) BoundaryConstraints = BoundaryConstraints()
%             end
%             if(nargin>1)
%                 noofval = numel(val);
%                 noofconst = numel(const);
%                 assert(noofval == noofconst || noofconst == 1, "Must be equal number of FitParameter objects and BoundaryConstraint objects");
% 
%                 values = val.multiGet("Value");
%                 bounds = const.multiGet("ScalarBounds");
% 
%                 [lb,ub] = bounds.multiGet("LowerBound","UpperBound");
%                 noofub = numel(ub);
%                 nooflb = numel(lb);
% 
%                 if(noofub > 0 && nooflb > 0)
%                     if(noofub > 1 && nooflb > 1)
%                         % array of lower bounds and upper bounds
%                         assert(noofub == nooflb,"Must have same number of lower and upper bounds")
%                         assert(noofub == noofval,"Each bound must correspond to a value")
%                     elseif(noofub > 1)
%                         % 1 lower bound, array of upper bounds
%                         assert(noofub == noofval,"Each bound must correspond to a value")
%                     else
%                         % 1 upper bound, array of lower bounds
%                         assert(nooflb == noofval,"Each bound must correspond to a value")
%                     end
% 
%                     checklower = all(lb<=values);
%                     checkupper = all(values<=ub);                
% 
%                     assert(checklower && checkupper,"Some Value is not between ScalarBounds")
%                 elseif (noofub > 0)
%                     % 1 or array of upper bounds, no lower bound
%                     assert(noofub == noofval || noofub == 1,"Each bound must correspond to a value")
%                     checkupper = all(values<=ub);                
% 
%                     assert(checkupper,"Some Value is not between ScalarBounds")
%                 elseif (nooflb > 0)
%                     % 1 or array of lower bounds, no upper bound
%                     assert(nooflb == noofval || nooflb == 1,"Each bound must correspond to a value")
%                     checklower = all(lb<=values);              
% 
%                     assert(checklower,"Some Value is not between ScalarBounds")
%                 end    
% 
%                 param(1:noofval) = NonLinearLeastSquaresParameter; 
%                 for i=1:noofval
%                     param(i).Value = val(i);
%                     if(noofconst==1)
%                         param(i).Constraints = const;
%                     else 
%                         param(i).Constraints = const(i);
%                     end    
%                 end 
%             end
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

