classdef BoundaryConstraints
    %BOUNDARYCONSTRAINTS Summary of this class goes here
    %   Detailed explanation goes here
    
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
    end
end
