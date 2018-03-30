classdef PassiveSpline1D < SplineInterpolant
    methods (Access = public)
        function Obj = PassiveSpline1D(varargin)
        % function Obj = PassiveSpline1D(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)
        % Class Constructor
            % Based on the sample points and the associated function values,
            % adjust the value of the terminal slopes.
            Obj.ROI_length = Obj.i_pts.getPtAt() - Obj.i_pts.getPtAt(1);
            Obj.MIN_EXTRAP_SLOPE = 0.05*max(abs(Obj.f_vals(Obj.colons{:}, :)), [], 1)/ROI_length;
            Obj.setupCoefficients()
        end
    end
end