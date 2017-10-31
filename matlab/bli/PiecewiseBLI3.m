classdef PiecewiseBLI3 < PiecewiseInterpolant
% Classdef PiecewiseBLI3
% Class definition for a 3D, Piecewise Barycentric Lagrange interpolant
% This reuses most of the code from the PiecewiseInterpolant class. The only
% significant difference is the declaration of the access handle BLI2
    methods (Access = protected)
        function populate(Obj)
            Obj.m_interp = squeeze(cell([Obj.n_pieces 1]));
            for interp_in = 1 : prod(Obj.n_pieces)
                Obj.m_interp{interp_in} = BLI3();
            end
        end
    end

    methods (Access = public)
        function Obj = PiecewiseBLI3(varargin)
        % function Obj = PiecewiseBLI3(f_handle, n_in_dims, bounds, order, i_type, smooth)
        % Class constructor
            Obj = Obj@PiecewiseInterpolant(varargin{:}, @BLI3);
            if ( isempty(varargin) )
                return;
            end
        end
    end
end
