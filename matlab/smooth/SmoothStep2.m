classdef SmoothStep2 < Interpolant
% SMOOTHSTEP2
% Extension of SmoothStep into 2 dimensions

    properties (Access = protected)
        i_type    = [];
        step_fit  = [];
        coeff_fit = [];
        % TODO: In order to use Save/Load functionality, we might have to change
        % the default value to SmoothStep().
    end

    methods (Access = public)
        function Obj = SmoothStep2(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % function Obj = SmoothStep2(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % Class constructor
            Obj = Obj@Interpolant(f_vals, in_dims, bounds, order, ...
                i_type_or_x_vals);

            % TODO: The line below (and the class member) is a quick hack to
            % get this working. There must be a better way of doing this.
            Obj.i_type = i_type_or_x_vals;
            Obj.fit([], []);
        end

        function coeffs = fit(Obj, vp_l, vp_r)
            % The FIT function in this case is a little tricky. In the 1D case,
            % we wanted to fit the derivative at 2 discrete points, now we have
            % to match the both partial derivatives (corresponding to the 2
            % dimensions) at a rectangular boundary

            % Express the function values as a SmoothStep(s) expression in the first dimension
            Obj.coeffs = zeros([Obj.op_dims, Obj.i_pts.n_pts(1), Obj.i_pts.n_pts(2)]);
            for p_i = 1:Obj.i_pts.n_pts(2)
                sstep_fit = BLI(Obj.f_vals(Obj.colons{:}, :, p_i), 1, ...
                    {Obj.bounds{1}}, Obj.order(1, 1), Obj.i_type);
                Obj.coeffs(Obj.colons{:}, :, p_i) = sstep_fit.coeffs;
            end

            Obj.step_fit = sstep_fit;
            Obj.coeff_fit = SmoothStep(Obj.coeffs, 1, {Obj.bounds{2}}, ...
                Obj.order(1, 2), Obj.i_type);
            % TODO: This is a recursive implementation which, in MATLAB, ends up
            % being extremely inefficient. Since the interpolant has a very
            % small number of coefficients, for speed, we can hardcode the
            % evaluation expression here instead of calling the underlying the
            % 1D interpolants recursively.
            if (nargout > 0)
                coeffs = Obj.coeffs;
            end
        end

        function [val, der] = computeWithDer(Obj, x_in)
            if (nargout > 1)
                der                          = zeros([Obj.op_dims, Obj.in_dims]);
                [coeffs, d_coeffs]           = Obj.coeff_fit.computeWithDer(x_in(2));
                [val, der(Obj.colons{:}, 1)] = Obj.step_fit.computeWithDer( ...
                    x_in(1), coeffs);
                der(Obj.colons{:}, 2)        = Obj.step_fit.computeWithDer( ...
                    x_in(1), d_coeffs);
                return;
            end
            l_coeffs = Obj.coeff_fit.computeWithDer(x_in(2));
            val      = Obj.step_fit.computeWithDer(x_in(1), l_coeffs);
        end

        function der = firstDerivativeAtPt(Obj, pt_index)
        
        end
function der = secondDerivativeAtPt(Obj, pt_index)

        end
    % end methods
    end
    
% end classdef
end
