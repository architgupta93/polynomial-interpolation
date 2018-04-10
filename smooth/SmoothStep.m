classdef SmoothStep < Interpolant
% SMOOTHSTEP
% This interpolant just creates a smooth patch between desired points

    properties (SetAccess = protected)
        
    end

    methods (Access = public)
        function Obj = SmoothStep(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % function Obj = SmoothStep(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % Class constructor
            Obj = Obj@Interpolant(f_vals, in_dims, bounds, order, ...
                i_type_or_x_vals);
            % TODO: Replace the argument list with varargin to support the use
            % of Save/Load functionality in Octave

            % The domain is [-1, 1]; rescaling the actual domain to this has
            % already been handled above.
            Obj.fit([], []);
        end 

        function coeffs = fit(Obj, vp_l, vp_r)
            % The scalar equations for this would be:
            %   p(x)   = ax3 + bx2 + cx + d
            %   p'(-1) = 3a - 2b + c        ... vp_l
            %   p'(1)  = 3a + 2b + c        ... vp_r
            %   p(-1)  = -a + b - c + d     ... v_l
            %   p(1)   =  a + b + c + d     ... v_r

            % Need to adjust the derivative for the scaling
            [~, dx]    = Obj.i_pts.rescaleShiftInput(0);
            if (isempty(vp_l))
                vp_l = zeros(Obj.op_dims);
            else
                vp_l = vp_l / dx;
            end

            if (isempty(vp_r))
                vp_r = zeros(Obj.op_dims);
            else
                vp_r = vp_r / dx;
            end

            np  = Obj.getNPts();
            v_l = Obj.f_vals(Obj.colons{:}, np);
            v_r = Obj.f_vals(Obj.colons{:}, 1);
            %   We have the following equations
            %   b      = (p'(1) - p'(-1))/4
            %   d      = -b + (p(1) + p(-1))/2
            %   a      = (p'(1) - p'(-1))/4 - (p(1) - p(-1))/4
            %   c      = -a + (p(1) - p(-1))/2
            % Create coefficients
            Obj.coeffs = zeros([Obj.op_dims 4]);

            Obj.coeffs(Obj.colons{:}, 2) = (vp_r - vp_l)/4.0;
            Obj.coeffs(Obj.colons{:}, 4) = - Obj.coeffs(Obj.colons{:}, 2) + ...
                (v_l + v_r)/2.0;
            Obj.coeffs(Obj.colons{:}, 1) =  (vp_r + vp_l)/4.0 - (v_r - v_l)/4.0;
            Obj.coeffs(Obj.colons{:}, 3) = - Obj.coeffs(Obj.colons{:}, 1) + ...
                (v_r - v_l)/2.0;

            if (nargout > 0)
                coeffs = Obj.coeffs;
            end
        end

        function [val, der] = computeWithDer(Obj, x_in, coeffs)
            if (nargin < 3)
                coeffs = Obj.coeffs;
            end
            n_pts          = Obj.getNPts();
            der            = zeros(Obj.op_dims);
            [x_in, dx_out] = Obj.i_pts.rescaleShiftInput(x_in);

            if (x_in < -1)
                val = Obj.f_vals(Obj.colons{:}, n_pts);
            elseif (x_in > 1)
                val = Obj.f_vals(Obj.colons{:}, 1);
            else
                % Compute the polynomial  terms
                x_in2    = times(x_in, x_in);
                x_in3    = times(x_in, x_in2);

                ext_dims = length(Obj.op_dims);
                x_vec    = shiftdim([x_in3, x_in2, x_in, 1], 1-ext_dims);
                val      = sum(x_vec .* coeffs, ext_dims+1);
            end
        end

        function der = firstDerivativeAtPt(Obj, pt_index)
            % TODO
            der = 0;
        end

        function der = secondDerivativeAtPt(Obj, pt_index)
            % TODO
            der = 0;
        end
    % end methods
    end
    
% end classdef
end
