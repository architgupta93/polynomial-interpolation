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
            Obj = Obj@PiecewiseInterpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            % Assigning a Chebyshev series object to handle each of the pieces
            % for interpolation

            varargin = {varargin{1:5}};   % Drop the last argument even it was supplied
            Obj.m_interp = cell(squeeze([Obj.n_pieces 1]));
            for p_i = 1:prod(Obj.n_pieces)
                l_bounds = Obj.getLocalBounds(p_i);
                varargin{3} = {l_bounds};
                %fprintf(2, ['DEBUG: Creating local BLI Interpolant with ', ...
                %    'bounds: \n']);
                %disp(l_bounds);
                Obj.m_interp{p_i} = BLI(varargin{:});
                % Obj.m_interp{p_i}.plotDerivatives(); %% DEBUG ONLY
                %if ( Obj.is_smooth && (p_i > 1) )
                %    % Remember that getPtAt(1) returns the last point for Chebyshev points as the points are ordered in
                %    % a decreasing order
                %    Obj.m_interp{p_i-1}.fit(Obj.m_interp{p_i}.firstDerivativeAtPt());
                %end
            end
        end

        function [val, der] = computeWithDer(Obj, x_in)
            if (nargout > 1)
                [val, der] = Obj.interpolate(x_in);
            else
                val = Obj.interpolate(x_in);
            end
        end
    end
end
