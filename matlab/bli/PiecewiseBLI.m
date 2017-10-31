classdef PiecewiseBLI < PiecewiseInterpolant
% PiecewiseBLI
% Class definition for a piecewise Barycentric-Lagrange interpolant
    methods (Access = protected)
        function populate(Obj)
            Obj.m_interp = squeeze(cell([Obj.n_pieces 1]));
            for interp_in = 1 : prod(Obj.n_pieces)
                Obj.m_interp{interp_in} = BLI();
            end
        end
    end

    methods (Access = public)
        function Obj = PiecewiseBLI(varargin)
        % function Obj = PiecewiseBLI(f_handle, n_in_dims, bounds, order, i_type, smooth)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(varargin{:}, @BLI);
        end
    end
end
