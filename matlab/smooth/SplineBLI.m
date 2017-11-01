classdef SplineBLI < BLI
% Classdef SplineBLI
% Class definition for smooth piecewise Barycentric-Lagrange interpolant
% Individual pieces (BLI) of interpolant are constructed by sampling a function
% in the appropriate domain and then they are stitched together using an
% overlaid spline of the interpolants.
    properties (Access = protected)
        spline_degree = 0;
        spline_coeff = 0;
    end

    methods (Access = public)
        function Obj = SplineBLI(f_vals, n_in_dims, bounds, order, i_type_or_x_vals, der)
            Obj = Obj@BLI(f_vals, n_in_dims, bounds, order, i_type_or_x_vals);
            
            if (nargin > 5)
                Obj.fit(der);
            end
        end

        function [val, der] = computeWithDer(Obj, x_in)
           [val, der] = computeWithDer@BLI(Obj, x_in);

           % Sorry for the code duplication
            x_vals = Obj.i_pts.getPts(1); 
            bpt = x_vals(end);
            x_rel = x_vals - x_in;  % x_in is a scalar and x_vals is a vector
            [found, loc] = builtin('_ismemberhelper', 0, x_rel);
            if (found)

            else
                val = val + ( Obj.spline_coeff * ( x_in - bpt) / (Obj.de_vals * (1./ x_rel)) );
            end
        end

        function fit(Obj, der)
            % Fitting the overlying spline.. In order to fit a smoothing term, we need to know the derivative at the
            % "RIGHT" boundary. It could very well have been the left derivative but we will stick to this convention
            % for the time being

            n_pts = Obj.getNPts(1);
            Obj.spline_degree = 1;
            piece_scale = Obj.getPtAt(1) - Obj.getPtAt(n_pts);
            Obj.spline_coeff = ( Obj.firstDerivativeAtPt(1) - der ) * ( Obj.wts(1) / piece_scale   );
        end

        function der = firstDerivative(Obj, x_in)
            der = 0; % TODO
        end

        function der = secondDerivative(Obj, x_in)
            der = 0; % TODO
        end
    end
end
