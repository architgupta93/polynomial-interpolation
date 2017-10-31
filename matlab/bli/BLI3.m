classdef BLI3 < Interpolant
% BLI3
% Implementation of a 3D Barycentric Lagrange Interpolant
    properties
        BLI_fit = BLI();
        nu_pols = BLI2();
    end

    methods
        function Obj = BLI3(varargin)
        % function Obj = BLI3(f_vals, in_dims, bounds, order, i_type_or_x_vals) 
        % Class constructor
            Obj = Obj@Interpolant(varargin{:});
            if (isempty(varargin))
                % For the Save/Load manual (octave) interface
                return;
            end

            % TODO: We should probably be using Obj.colons here somewhere
            g_nu = zeros(Obj.op_dims, Obj.i_pts.n_pts(1), ...
                Obj.i_pts.n_pts(2), Obj.i_pts.n_pts(3));
            
            for p_i = 1:Obj.i_pts.n_pts(3)
                for p_j = 1:Obj.i_pts.n_pts(2)
                    BLI_fit = BLI(Obj.f_vals(Obj.colons{:}, :, p_j, p_i), 1, ...
                        {Obj.bounds{1}}, Obj.order(1, 1), ...
                        varargin{end});
                    %   ^^ This is i_type_or_x_vals

                    g_nu(:, :, p_j, p_i) = BLI_fit.nu_vals;
                end
            end

            Obj.BLI_fit = BLI_fit;
            Obj.nu_pols = BLI2(g_nu, 2, {Obj.bounds{2}, Obj.bounds{3}}, ...
                Obj.order(1, 2:3), varargin{end});
        end
        
        % API Function(s) for evaluation
        function [val, der] = computeWithDer(Obj, x_in)
            if (nargout > 1)
                der = zeros([Obj.op_dims, Obj.in_dims]);
                [i_nuvals, d_nuvals] = Obj.nu_pols.computeWithDer(x_in(2:3));
                [val, der(Obj.colons{:}, 1)]  = Obj.BLI_fit.computeWithDer(x_in(1), i_nuvals);
                der(Obj.colons{:}, 2)  = Obj.BLI_fit.computeWithDer(x_in(1), ...
                    d_nuvals(Obj.colons{:}, :, 1));
                der(Obj.colons{:}, 3)  = Obj.BLI_fit.computeWithDer(x_in(1), ...
                    d_nuvals(Obj.colons{:}, :, 2));
                return;
            end

            val = Obj.BLI_fit.computeWithDer(x_in(1), ...
                Obj.nu_pols.computeWithDer(x_in(2:3)));
        end

        % Some other rarely used access functions
        function der = secondDerivative(Obj, x_in)
            der = 0;
            % TODO
        end

        function der = secondDerivativeAtPt(Obj, pc_index)
            der = 0;
            % TODO
        end

        function der = firstDerivativeAtPt(Obj, pc_index)
            der = 0;
            % TODO
        end
    % end methods
    end
    
% end classdef
end
