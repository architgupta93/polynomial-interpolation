classdef PassiveSpline1D < SplineInterpolant

    methods (Access = protected)
        function [b, M] = adjustForBoundarySlope(Obj, h, b, M)
            mid_point_idx   = ceil(Obj.i_pts.getNPts()/2);
            mid_point_val   = Obj.f_vals(Obj.colons{:}, mid_point_idx);
            MIN_SLOPE_RIGHT = -(Obj.f_vals(Obj.colons{:}, end) - mid_point_val).*Obj.MIN_EXTRAP_SLOPE;
            MIN_SLOPE_LEFT  = (Obj.f_vals(Obj.colons{:}, 1) -  mid_point_val).*Obj.MIN_EXTRAP_SLOPE;
            
                                % Minimum value of the absolute slope that our
                                % interpolation should have at the boundary.
                                % This ensures that Newton lands you in the
                                % right area in a single step
            
            % Adjusting b so that the interpolation has the correct slope (and
            % the minimum value) on the right/left hand side

            n_plus_1 = Obj.i_pts.getNPts(1);
            b(:, end) = b(:, end) + 3*(Obj.f_vals(Obj.colons{:}, end)-Obj.f_vals(Obj.colons{:}, end-1))/h(1,end) ...
                - 3*MIN_SLOPE_RIGHT;  % Notice that this automatically gives us
                                        % the correct sign for f'(x)

            M(n_plus_1-2,n_plus_1-2) = M(n_plus_1-2,n_plus_1-2) - h(1,end)/2;

            b(:,1) = b(:,1) + 3*(Obj.f_vals(Obj.colons{:}, 1)-Obj.f_vals(Obj.colons{:}, 2))/h(1,1) ...
                + 3*MIN_SLOPE_LEFT;

            M(1,1) = M(1,1) - h(1,1)/2;
        end
    end

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
            % Obj.MIN_EXTRAP_SLOPE = 0.5*max(abs(Obj.f_vals(Obj.colons{:}, :)), [], 1+length(Obj.colons))/ROI_length;
            Obj.MIN_EXTRAP_SLOPE = 2.0/ROI_length;
            Obj.setupCoefficients()
        end
    end
end