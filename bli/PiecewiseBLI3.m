classdef PiecewiseBLI3 < PiecewiseInterpolant
% Classdef PiecewiseBLI3
% Class definition for a 3D, Piecewise Barycentric Lagrange interpolant
% This reuses most of the code from the PiecewiseInterpolant class. The only
% significant difference is the declaration of the access handle BLI2
    methods (Access = public)
        function Obj = PiecewiseBLI3(f_handle, n_in_dims, bounds, order, i_type)
        % function Obj = PiecewiseBLI3(f_handle, n_in_dims, bounds, order, i_type)
        % Class constructor
            % TODO: Check that n_in_dims is actually sane. 
            % It was observed recently that you could take a 3 input function
            % and pass it into PiecewiseBLI2 and it would actually produce some
            % interpolant. This is very dangerous!
            Obj = Obj@PiecewiseInterpolant(f_handle, n_in_dims, bounds, order, i_type, @BLI3);
        end
    end
end
