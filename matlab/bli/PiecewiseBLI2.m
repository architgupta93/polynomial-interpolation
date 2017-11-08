classdef PiecewiseBLI2 < PiecewiseInterpolant
% Classdef PiecewiseBLI2
% Class definition for a 2D, Piecewise Barycentric Lagrange interpolant
% This reuses most of the code from the PiecewiseInterpolant class. The only
% significant difference is the declaration of the access handle BLI2
    methods (Access = public)
        function Obj = PiecewiseBLI2(varargin)
        % function Obj = PiecewiseBLI2(f_handle, n_in_dims, bounds, order, i_type)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(varargin{:}, @BLI2);
        end
    end
end
