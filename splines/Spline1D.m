classdef Spline1D < SplineInterpolant
    methods (Access = public)
        function Obj = Spline1D(varargin)
        % function Obj = Spline1D(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)        % Class Constructor
            Obj = Obj@SplineInterpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            Obj.setupCoefficients();
        end
    end
end
