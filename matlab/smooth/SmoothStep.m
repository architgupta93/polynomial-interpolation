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

            % Create coefficients
            Obj.coeffs = zeros([Obj.op_dims 2]);
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
                % Let's start with an actual step function for which the value
                % of the derivative at the boundary points is 0.
                %
                % The domain is [-1, 1]; rescaling the actual domain to this has
                % already been handled above.
                %
                % The scalar equations for this would be:
                %   p(x)   = ax3 + bx2 + cx + d
                %   p'(-1) = 3a - 2b + c
                %   p'(1)  = 3a + 2b + c
                %   p(-1)  = -a + b - c + d 
                %   p(1)   =  a + b + c + d
                %
                % Zero derivatives at the boundaries already gives us:
                %   b      = 0
                %   c      = -3a
                %   d      = (p(1) + p(-1))/2
                %   a      = (d - p(1))/2
            end
        end
    % end methods
    end
    
% end classdef
end
