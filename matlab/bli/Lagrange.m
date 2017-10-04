classdef Lagrange < Interpolant
    properties (SetAccess = protected, GetAccess = public)
        pols = [];
    end

    methods (Access = public)
        function Obj =  Lagrange(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)

            Obj = Obj@Interpolant(f_vals, n_in_dims, bounds, order, i_type_or_x_vals);

            colons = cell( size(Obj.op_dims) );
            [colons{:}] = deal(':');
            % Initializing all the lagrange polynomials
            wts  = zeros([Obj.op_dims, Obj.i_pts.getNPts(1)   ]);
            Obj.pols = zeros([Obj.op_dims, Obj.i_pts.getNPts(1)   ]);
            pols = zeros(Obj.i_pts.getNPts(1));
            all_pts = Obj.i_pts.pts{1};
            diff_dims = ndims(pols) - ndims(wts);
            for p_i = 1 : Obj.i_pts.getNPts(1)
                pols(p_i, :) = poly( [all_pts(1:p_i-1); all_pts(p_i+1:end)] );
                wts(colons{:}, p_i) = rdivide(Obj.f_vals(colons{:}, p_i), polyval(pols(p_i, :), Obj.i_pts.pts{1}(p_i)));
                Obj.pols = Obj.pols + times( shiftdim(pols(p_i, :), diff_dims), wts(colons{:}, p_i) );
            end
        end

        function [val, der] = computeWithDer(Obj, x_in, coeffs)
            [x_in, dx_out] = Obj.i_pts.rescaleShiftInput(x_in);
            if ( nargin < 3 )
                coeffs = Obj.pols;
            end

            if (isscalar(Obj.op_dims))
                val = zeros([Obj.op_dims 1]);
                der = zeros([Obj.op_dims 1]);
            else
                val = zeros(Obj.op_dims);
                der = zeros(Obj.op_dims);
            end
            sub_i = cell( size(Obj.op_dims) );
            for p_i = 1 : prod(Obj.op_dims)
                [sub_i{:}] = ind2sub( Obj.op_dims, p_i );
                val(sub_i{:}) = polyval( squeeze(coeffs(sub_i{:}, :)), x_in );
                der(sub_i{:}) = dx_out * polyval( polyder(squeeze(coeffs(sub_i{:}, :))), x_in );
            end
        end

        function der = firstDerivativeAtPt(Obj, pt_index)
            % This function can be used to obtain the value of the first derivative at one of sample points. There is a
            % separate function for this as the formula can be somewhat easier than evaluating the derivative at an
            % arbitrary point within the sampling domain

            if ( nargin < 2)
                pt_index = Obj.getNPts(1);
            end

            [ ~, der] = Obj.computeWithDer( Obj.i_pts.getPtAt(pt_index) );
        end

        function der = secondDerivativeAtPt(Obj, pt_index)
            der = 0;
            % TODO: Not sure if this will be required
        end

        function der = secondDerivative(Obj, x_in)
            der = 0;
            % TODO: Again not sure if it will ever been needed but might be good to have it
        end
    end
end
