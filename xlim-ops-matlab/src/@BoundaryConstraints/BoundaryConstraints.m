classdef BoundaryConstraints
    % BOUNDARYCONSTRAINTS Holds boundary constraints to be used with a particular parameter
    
    properties (Dependent)
        % SCALARBOUNDS - Lower/upper bounds around parameter
        ScalarBounds (1,1) ScalarParameterBounds
        % LINEARINEQUALITYVECTOR -  Row vector of bounds to apply for linear inequality constraints
        LinearInequalityVector (1,:) {mustBeNumeric}
        % LINEARINEQUALITYSCALAR - Scalar bound to apply for linear inequalities
        LinearInequalityScalar (1,1) {mustBeNumeric}
        % LINEAREQUALITYVECTOR - Row vector of bounds to apply for linear equality constraints
        LinearEqualityVector (1,:) {mustBeNumeric}
        % LINEAREQUALITYSCALAR - Scalar bound to apply for linear equalities 
        LinearEqualityScalar (1,1) {mustBeNumeric}
        % NONLINEARCONSTRAINTS - Function handle to apply for non-linear constraints
        NonLinearConstraints (1,1) function_handle
    end

    properties (Access = protected,Hidden)
        % HIDDENSCALARBOUNDS
        HiddenScalarBounds
        % HIDDENLINEARINEQUALITYVECTOR
        HiddenLinearInequalityVector
        % HIDDENLINEARINEQUALITYSCALAR
        HiddenLinearInequalityScalar
        % HIDDENLINEAREQUALITYVECTOR
        HiddenLinearEqualityVector
        % HIDDENLINEAREQUALITYSCALAR
        HiddenLinearEqualityScalar
        % HIDDENNONLINEARCONSTRAINTS
        HiddenNonLinearConstraints
    end
    
    methods
        function const = BoundaryConstraints(args)
            %BOUNDARYCONSTRAINTS - Holds boundary constraints to be used with a particular parameter
            arguments
                args.?BoundaryConstraints
            end
            arc = namedargs2cell(args);
            for i = 1:2:numel(arc)
                const.(arc{i}) = arc{i + 1};
            end
        end
        
        function const = set.ScalarBounds(const,scalarbounds)
            const.HiddenScalarBounds = scalarbounds;
        end

        function scalarbounds = get.ScalarBounds(const)
            scalarbounds = const.HiddenScalarBounds;
        end

        function const = set.LinearInequalityVector(const,vector)
            const.HiddenLinearInequalityVector = vector;
        end

        function vector = get.LinearInequalityVector(const)
            vector = const.HiddenLinearInequalityVector;
        end

        function const = set.LinearInequalityScalar(const,scalar)
            const.HiddenLinearInequalityScalar = scalar;
        end

        function scalar = get.LinearInequalityScalar(const)
            scalar = const.HiddenLinearInequalityScalar;
        end

        function const = set.LinearEqualityVector(const,vector)
            const.HiddenLinearEqualityVector = vector;
        end

        function vector = get.LinearEqualityVector(const)
            vector = const.HiddenLinearEqualityVector;
        end

        function const = set.LinearEqualityScalar(const,scalar)
            const.HiddenLinearEqualityScalar = scalar;
        end

        function scalar = get.LinearEqualityScalar(const)
            scalar = const.HiddenLinearEqualityScalar;
        end

        function const = set.NonLinearConstraints(const,nonlinconst)
            assert(nargin(nonlinconst) == 1);
            assert(nargout(nonlinconst) == 2);
            const.HiddenNonLinearConstraints = nonlinconst;
        end

        function nonlinconst = get.NonLinearConstraints(const)
            nonlinconst = const.HiddenNonLinearCOnstraints;
        end

        function varargout = multiGet(multibound, prop)
            %Property parser for single or array of BoundaryConstraints objects
            arguments
                multibound (:,1) BoundaryConstraints
            end

            arguments (Repeating)
                prop (1,1) {mustBeMember(prop,["ScalarBounds","LinearInequalityVector","LinearInequalityScalar","LinearEqualityVector","LinearEqualityScalar","NonLinearConstraints"])}
            end

            varargout = cell(size(prop));

            for i=1:numel(prop)
                switch prop{i}
                    case {"LinearInequalityVector","LinearEqualityVector"} 
                        try 
                            varargout{i} = vertcat(multibound.(prop{i}));
                        catch 
                            error("Some Linear Vectors are not of the same length")
                        end 
                    otherwise
                        varargout{i} = vertcat(multibound.(prop{i}));
                end
            end
        end 

        function multiconst = multiSet(multiconst,args)
            arguments
                multiconst (:,1) BoundaryConstraints
                args.ScalarBounds (:,1) ScalarParameterBounds
                args.LinearInequalityVector (:,:) {mustBeNumeric}
                args.LinearInequalityScalar (:,1) {mustBeNumeric}
                args.LinearEqualityVector (:,:) {mustBeNumeric}
                args.LinearEqualityScalar (:,1) {mustBeNumeric}
                args.NonLinearConstraints (:,1) function_handle
            end
            
            nbounds = numel(multiconst);

            isscalarb = isfield(args,"ScalarBounds");
            if(isscalarb)
                nscalarb = numel(args.ScalarBounds);
                isscalarb1 = nscalarb == 1;
                assert(nscalarb == nbounds || isscalarb1, "No of ScalarBounds and no of BoundaryConstraint objects do not match")
            end

            isLinInVec = isfield(args,"LinearInequalityVector");
            if(isLinInVec)
                nlininvec = size(args.LinearInequalityVector,1);
                islininvec1 = nlininvec == 1;
                assert(nlininvec == nbounds || islininvec1, "No of LinearInequalityVectors and no of BoundaryConstraint objects do not match")
            end

            isLinEqVec = isfield(args,"LinearEqualityVector");
            if(isLinEqVec)
                nlineqvec = size(args.LinearEqualityVector,1);
                islineqvec1 = nlineqvec == 1;
                assert(nlineqvec == nbounds || islineqvec1, "No of LinearEqualityVectors and no of BoundaryConstraint objects do not match")
            end

            isLinInScal = isfield(args,"LinearInequalityScalar");
            if(isLinInScal)
                nlininscal = numel(args.LinearInequalityScalar);
                islininscal1 = nlininscal == 1;
                assert(nlininscal == nbounds || islininscal1, "No of LinearInequalityScalars and no of BoundaryConstraint objects do not match")
            end

            isLinEqScal = isfield(args,"LinearEqualityScalar");
            if(isLinEqScal)
                nlineqscal = numel(args.LinearEqualityScalar);
                islineqscal1 = nlineqscal == 1;
                assert(nlineqscal == nbounds || islineqscal1, "No of LinearEqualityScalars and no of BoundaryConstraint objects do not match")
            end

            isnonlin = isfield(args,"NonLinearConstraints");
            if(isnonlin)
                nnonlin = numel(args.NonLinearConstraints);
                isnonlin1 = nnonlin == 1;
                assert(nnonlin == nbounds || isnonlin1, "No of NonLinearConstraints provided and no of BoundaryConstraint objects do not match")
            end

            for i = 1:nbounds
                if(isscalarb)
                    if (~isscalarb1)
                        multiconst(i).HiddenScalarBounds = args.ScalarBounds(i);
                    else
                        multiconst(i).HiddenScalarBounds = args.ScalarBounds;
                    end
                end
                if(isLinInVec)
                    if (~islininvec1)
                        multiconst(i).HiddenLinearInequalityVector = args.LinearInequalityVector(i,:);
                    else
                        multiconst(i).HiddenLinearInequalityVector = args.LinearInequalityVector;
                    end
                end
                if(isLinEqVec)
                    if (~islineqvec1)
                        multiconst(i).HiddenLinearEqualityVector = args.LinearEqualityVector(i,:);
                    else
                        multiconst(i).HiddenLinearEqualityVector = args.LinearEqualityVector;
                    end
                end
                if(isLinInScal)
                    if (~islininscal1)
                        multiconst(i).HiddenLinearInequalityScalar = args.LinearInequalityScalar(i);
                    else
                        multiconst(i).HiddenLinearInequalityScalar = args.LinearInequalityScalar;
                    end
                end
                if(isLinEqScal)
                    if (~islineqscal1)
                        multiconst(i).HiddenLinearEqualityScalar = args.LinearEqualityScalar(i);
                    else
                        multiconst(i).HiddenLinearEqualityScalar = args.LinearEqualityScalar;
                    end
                end
                if(isnonlin)
                    if (~isnonlin1)
                        multiconst(i).HiddenNonLinearConstraints = args.NonLinearConstraints(i);
                    else
                        multiconst(i).HiddenNonLinearConstraints = args.NonLinearConstraints;
                    end
                end
            end 
        end
    end
end

