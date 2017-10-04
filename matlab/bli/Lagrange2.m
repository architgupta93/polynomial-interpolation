classdef Lagrange2 < Interpolant
    properties (SetAccess = protected, GetAccess = public)
        pols = struct();
        g_pols = [];
        Lg_fit = struct();
    end

    methods (Access = public)
        function Obj = Lagrange2(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)
            Obj = Obj@Interpolant(f_vals, n_in_dims, bounds, order, i_type_or_x_vals);

            % Express the function values as a Lagrange interpolant in the first dimension
            Obj.g_pols = zeros(Obj.op_dims, Obj.getNPts(1), Obj.getNPts(2));
            for p_i = 1:Obj.getNPts(2)
                Lg_fit = Lagrange( squeeze(Obj.f_vals(:, :, p_i)), 1, bounds(:,1), order(1, 1), i_type_or_x_vals);
                Obj.g_pols(:, :, p_i) = Lg_fit.pols;
            end
            
            % In the last step, we took all the 'x' values and expressed the function as a Lagrange Interpolant for each
            % value of 'x'. Now, we need to take the polynomial coefficients for each of these interpolants and express
            % them as a Lagrange Interpolant in the other dimension.
            Obj.Lg_fit = Lg_fit;
            Obj.pols = Lagrange(Obj.g_pols, 1, bounds(:,2), order(1, 2), i_type_or_x_vals);
        end

        function der = firstDerivativeAtPt(Obj, pc_index)
            % TODO
            der = zeros(Obj.n_in_dims, 1);
        end

        function der = secondDerivativeAtPt(Obj, pc_index)
            % TODO
            der = zeros(Obj.n_in_dims, 1);
        end

        function [val, der] = computeWithDer(Obj, x_in)
            [i_fvals, i_dfvals] = Obj.pols.computeWithDer(x_in(2));
            der = zeros([Obj.op_dims, Obj.n_in_dims]);
            [val, der(Obj.colons{:}, 1)] = Obj.Lg_fit.computeWithDer(x_in(1), i_fvals);
            der(Obj.colons{:}, 2) = Obj.Lg_fit.computeWithDer(x_in(1), i_dfvals);
        end

        function der = secondDerivative(Obj, x_in)
            der = 0; % TODO
        end

        function y_out = evaluate(Obj, x_in)
            i_fvals = Obj.pols.evaluate(x_in(2));
            y_out = Obj.Lg_fit.evaluateWithCoeffs(x_in(1), i_fvals);
        end
    end
end
