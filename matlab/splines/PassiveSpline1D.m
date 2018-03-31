classdef PassiveSpline1D < SplineInterpolant
    methods (Access = public)
        function Obj = PassiveSpline1D(varargin)
        % function Obj = PassiveSpline1D(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)
        % Class Constructor
            Obj = Obj@SplineInterpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            % Based on the sample points and the associated function values,
            % adjust the value of the terminal slopes.
            ROI_length = Obj.i_pts.getPtAt(-1) - Obj.i_pts.getPtAt(1);
            Obj.MIN_EXTRAP_SLOPE = 0.1*max(abs(Obj.f_vals(Obj.colons{:}, :)), [], 1+length(Obj.colons))/ROI_length;
            Obj.setupCoefficients()
        end
    end
end