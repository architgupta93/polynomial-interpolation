classdef PiecewiseSpline < PiecewiseInterpolant
% PiecewiseSpline
% Class definition for a piecewise <Spline> interpolant
    methods (Access = public)
        function Obj = PiecewiseSpline(f_handle, n_in_dims, bounds, order, i_type, smooth)
        % function Obj = PiecewiseSpline(f_handle, n_in_dims, bounds, order, i_type, smooth)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(f_handle, n_in_dims, bounds, order, i_type, @Spline1D);
        end
    end
end
