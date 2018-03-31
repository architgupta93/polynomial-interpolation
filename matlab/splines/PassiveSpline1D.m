classdef PassiveSpline1D < SplineInterpolant

    methods (Access = protected)
        function [z, M] = adjustForBoundarySlope(Obj, h, b, M)
            % mid_point_idx   = ceil(Obj.i_pts.getNPts()/2);
            % mid_point_val   = Obj.f_vals(Obj.colons{:}, mid_point_idx);

            % Keep in mind that our sample points (and thus sample values)
            % are produced in reverse. The sample point corresponding to the
            % last index is actually the first value.
            n_plus_1 = Obj.i_pts.getNPts(1);
            % MIN_SLOPE_RIGHT = sign(Obj.f_vals(Obj.colons{:}, n_plus_1)).*Obj.MIN_EXTRAP_SLOPE;
            % MIN_SLOPE_LEFT  = -sign(Obj.f_vals(Obj.colons{:}, 1)).*Obj.MIN_EXTRAP_SLOPE;
            MIN_SLOPE_LEFT  = 0;
            MIN_SLOPE_RIGHT = 0;
            
                                % Minimum value of the absolute slope that our
                                % interpolation should have at the boundary.
                                % This ensures that Newton lands you in the
                                % right area in a single step
            
            % Adjusting b, M so that the interpolation has the correct slope (and
            % the minimum value) on the right/left hand side
            b(:, end) = b(:, end) + 3*(Obj.f_vals(Obj.colons{:}, n_plus_1)-Obj.f_vals(Obj.colons{:}, n_plus_1-1))/h(1,end) ...
                - 3*MIN_SLOPE_RIGHT;
            b(:,1) = b(:,1) + 3*(Obj.f_vals(Obj.colons{:}, 1)-Obj.f_vals(Obj.colons{:}, 2))/h(1,1) ...
                + 3*MIN_SLOPE_LEFT;
            M(n_plus_1-2,n_plus_1-2) = M(n_plus_1-2,n_plus_1-2) - h(1,end)/2;
            M(1,1) = M(1,1) - h(1,1)/2;

            M_inv = full( inv(M) );
            z = zeros( size(b) );
            diff_dims = ndims(h) - ndims(Obj.f_vals);
            for c_index = 1 : size(M, 2)
                z(Obj.colons{:}, c_index) = sum( times(b, shiftdim(M_inv(:, c_index)', diff_dims)), ndims(b) );
            end

            z1 = z(Obj.colons{:}, 1);
            zn_minus_1 = z(Obj.colons{:}, n_plus_1-2);

            s_z  = size(z);
            nd_z = length(s_z);
            z = cat(nd_z, 3*(Obj.f_vals(Obj.colons{:}, 2)-Obj.f_vals(Obj.colons{:}, 1))/h(1,1)^2 - z1/2.0 - 3.0*MIN_SLOPE_LEFT/h(1,1), ...
                z, ...
                3*(Obj.f_vals(Obj.colons{:}, n_plus_1-1)-Obj.f_vals(Obj.colons{:}, n_plus_1))/h(end,1)^2 - zn_minus_1/2.0 + 3.0*MIN_SLOPE_RIGHT/h(end,1)); 
                % Natural spline assumes that the second derivative at the two
                % end points is 0, which has been put in here
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
            ROI_length = Obj.i_pts.getPtAt(1) - Obj.i_pts.getPtAt(-1);
            % Obj.MIN_EXTRAP_SLOPE = 0.5*max(abs(Obj.f_vals(Obj.colons{:}, :)), [], 1+length(Obj.colons))/ROI_length;
            Obj.MIN_EXTRAP_SLOPE = 2.0/ROI_length;
            Obj.setupCoefficients()
        end
    end
end