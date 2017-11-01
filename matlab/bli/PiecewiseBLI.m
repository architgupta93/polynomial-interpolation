classdef PiecewiseBLI < PiecewiseInterpolant
% PiecewiseBLI
% Class definition for a piecewise Barycentric-Lagrange interpolant
    methods (Access = public)
        function Obj = PiecewiseBLI(f_handle, n_in_dims, bounds, order, ...
            i_type, smooth)
        % function Obj = PiecewiseBLI(f_handle, n_in_dims, bounds, order,
        %     i_type, smooth)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(f_handle, n_in_dims, bounds, order, i_type);
            if (nargin < 6)
                smooth = false;
            end

            Obj.is_smooth = smooth;
            if (smooth)
                Obj.setAccessHandle(@SplineBLI, f_handle);
            else
                Obj.setAccessHandle(@BLI, f_handle);
            end
        end
    end
end
