classdef ChebyshevPoints2Analyst < ChebyshevPoints2
    methods (Access = public)
        function obj = ChebyshevPoints2Analyst(n_in_dims, n_pts, bounds)
            obj = obj@ChebyshevPoints2(n_in_dims, n_pts, bounds);    
        end

        function coeff_vals = getSeriesCoeffs(Obj, f_vals)
        % Express the passed function values f_vals as a chebyshev series.
        % Internally, we are trying to represent the entire function f (not its
        % discretized values) as a chebyshev series. However, this can be done
        % efficiently using DCT if one is aware of the discrete orthogonality of
        % the Chebyshev Points. TODO: Derive the discrete orthogonality for
        % Chebyshev points of type 2
            if (Obj.roots_or_extremas) % Using roots
                coeff_vals = DCT1(f_vals, 'CHEB_SERIES');
            else % Use extremas
                coeff_vals = DCT2(f_vals, 'CHEB_SERIES');
            end
        end
    end
end
