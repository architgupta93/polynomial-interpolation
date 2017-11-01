classdef PiecewiseBLI < PiecewiseInterpolant
% PiecewiseBLI
% Class definition for a piecewise Barycentric-Lagrange interpolant
    methods (Access = public)
        function Obj = PiecewiseBLI(f_handle, n_in_dims, bounds, order, i_type)
        % function Obj = PiecewiseBLI(f_handle, n_in_dims, bounds, order, i_type)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(f_handle, n_in_dims, bounds, order, i_type, @BLI);
        end
    end
end
