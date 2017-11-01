classdef PiecewiseSpline < PiecewiseInterpolant
% PiecewiseSpline
% Class definition for a piecewise <Spline> interpolant
    methods (Access = protected)
        function populate(Obj)
            Obj.m_interp = squeeze(cell([Obj.n_pieces 1]));
            for interp_in = 1 : prod(Obj.n_pieces)
                Obj.m_interp{interp_in} = Spline1D();
            end
        end
    end

    methods (Access = public)
        function Obj = PiecewiseSpline(varargin)
        % function Obj = PiecewiseSpline(f_handle, n_in_dims, bounds, order, i_type, smooth)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(varargin{:}, @Spline1D);
        end
    end
end
