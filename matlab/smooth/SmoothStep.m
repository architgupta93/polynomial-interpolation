classdef SmoothStep < Interpolant
% SMOOTHSTEP
% This interpolant just creates a smooth patch between desired points

    properties (SetAccess = protected)
        
    end

    methods (Access = public)
        function Obj = SmoothStep(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % function Obj = Smooth(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % Class constructor
            Obj = Obj@Interpolant(f_vals, in_dims, bounds, order, ...
                i_type_or_x_vals);

            % The domain is [-1, 1]; rescaling the actual domain to this has
            % already been handled above.
            Obj.generateCoeffs([], []);
        end 

        function coeffs = generateCoeffs(Obj, vp_l, vp_r)
            % The scalar equations for this would be:
            %   p(x)   = ax3 + bx2 + cx + d
            %   p'(-1) = 3a - 2b + c        ... vp_l
            %   p'(1)  = 3a + 2b + c        ... vp_r
            %   p(-1)  = -a + b - c + d     ... v_l
            %   p(1)   =  a + b + c + d     ... v_r
            if (isempty(vp_l))
                vp_l = zeros(Obj.op_dims);
            end

            if (isempty(vp_r))
                vp_r = zeros(Obj.op_dims);
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
            Obj.coeffs(Obj.colons{:}, 4) = -Obj.coeffs(Obj.colons{:}, 2) + (v_l + v_r)/2.0;
            Obj.coeffs(Obj.colons{:}, 1) =  Obj.coeffs(Obj.colons{:}, 2) - (v_l - v_r)/4.0;
            Obj.coeffs(Obj.colons{:}, 3) = -Obj.coeffs(Obj.colons{:}, 1) + ...
                (v_r - v_l)/2.0;

            if (nargout > 0)
                coeffs = Obj.coeffs;
            end
        end

        function [val, der] = computeWithDer(Obj, x_in)
            n_pts          = Obj.getNPts();
            der            = zeros(Obj.op_dims);
            [x_in, dx_out] = Obj.i_pts.rescaleShiftInput(x_in);

            if (x_in < Obj.getPtAt(n_pts))
                val = Obj.f_vals(Obj.colons{:}, n_pts);
            elseif (x_in < Obj.getPtAt(n_pts))
                val = Obj.f_vals(Obj.colons{:}, 1);
            else
                % Compute the polynomial  terms
                x_in2    = times(x_in, x_in);
                x_in3    = times(x_in, x_in2);


                ext_dims = length(Obj.op_dims);
                x_vec    = shiftdim([x_in3, x_in2, x_in, 1], 1-ext_dims);
                val      = dx_out * sum(x_vec .* ...
                            Obj.coeffs, ext_dims+1);
            end
        end
    % end methods
    end
    
% end classdef
end
