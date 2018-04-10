classdef ChebyshevPoints1Analyst < ChebyshevPoints1
    methods (Access = protected, Static = true)
        function p_val = evaluate1DPolynomial(x_in, degree)
            p_val = zeros(degree, 1);
            p_val(1) = 1;
            p_val(2) = x_in;

            i = 2;
            while ( i < degree ) % Polynomial degree
                p_val(i+1) = 2 * ( x_in .* p_val(i)) - p_val(i-1);
                i = i+1;
            end
        end
    end

    methods (Access = public)
        function Obj = ChebyshevPoints1Analyst(n_in_dims, order, bounds);
            Obj = Obj@ChebyshevPoints1(n_in_dims, order, bounds);
        end

        function coeff_vals = getSeriesCoeffs(Obj, f_vals)
        % Express the passed function values f_vals as a chebyshev series.
        % Internally, we are trying to represent the entire function f (not its
        % discretized values) as a chebyshev series. However, this can be done
        % efficiently using DCT if one is aware of the discrete orthogonality of
        % the Chebyshev Points. TODO: Derive the discrete orthogonality for
        % Chebyshev points of type 2

            % Formulation: Let us call the polynomials T_i(x). Because of the
            % discrete orthogonality relation, we get:
            %
            %       sum''( T_r(x_j) * T_s(x_j) ) = K_r & delta_rs;
            %
            % Look at "Numerical methods for special functions, page 59", or
            % Notes on Chebyshev expansions for the meaning of various
            % quantities here. If we can express f as a truncated Chebyshev
            % series of the first 'n' Chebyshev polynomials, we get:
            %
            %                   f(x) = sum( c_r * T_r(x) )
            %
            % Using the discrete orthogonality at the extremas of T_n(x), we get
            %
            %       sum''( f(x_j) * T_s(x_j) ) = sum''( T_r(x_j) * T_s(x_j) ) 
            %
            % The exact transform to be taken depends on which sample points are
            % being used. For zeros of T_(n+1), we should use DCT2. Whereas, for
            % extremas of T_n, we should use DCT1

            % UPDATE: October 4, 2017
            % Because of the changes in which f_vals are stored inside
            % interpolants, these are noew row vectors, i.e., the last index
            % corresponds to the sample index. For this, we need to take a
            % transpose before passign on for DCT
            if (Obj.roots_or_extremas)  % Using roots
                coeff_vals = amateur__DCT2(f_vals, 'CHEB_SERIES');
            else    % Using extremas
                coeff_vals = amateur__DCT1(f_vals, 'CHEB_SERIES');
            end
        end

        function p_val = evaluatePolynomial(Obj, x_in)
            % Before we do any computation with x_in, we need to convert it to
            % the local coordinates and also scale it using the current center.
            % For example, if the current value for "bounds" in dimension i is
            % [0; 1], then a  value of x_in(i) = 0.5 actually translates to a
            % value 0 for chebyshev polynomials

            x_in = Obj.rescaleShiftInput(x_in);

            % For each element in x_in, we will compute a vector of evaluated
            % polynomials and return this inside a cell array to lower the
            % communication overhead. The tensor product of these vectors,
            % followed by scaling with coefficients and summing can be done in
            % the class that has called this function
            p_val = cell(Obj.n_in_dims, 1);
            for d_i = 1:Obj.n_in_dims
                p_val{d_i, 1} = ChebyshevPoints1Analyst.evaluate1DPolynomial(...
                    x_in(d_i), Obj.n_pts(d_i));
            end
        end
    end
end
