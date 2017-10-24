classdef BLI2 < Interpolant
    properties (Access = protected)
        BLI_fit = BLI();
        nu_pols = BLI();
    end

    methods (Access = public)
        function Obj = BLI2(varargin)
        % function Obj = BLI2(f_vals, in_dims, bounds, order, i_type_or_x_vals) 
            Obj = Obj@Interpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            % Express the function values as a BLI expression in the first dimensioN
            BLI_fit = cell(Obj.i_pts.n_pts(2), 1);
            g_nu = zeros(Obj.op_dims, Obj.i_pts.n_pts(1), Obj.i_pts.n_pts(2));
            
            for p_i = 1:Obj.i_pts.n_pts(2)
                BLI_fit = BLI( squeeze(Obj.f_vals(:, :, p_i)), 1, Obj.bounds{1}, Obj.order(1, 1), varargin{end});
                %                                                           This is i_type_or_x_vals ^^
                g_nu(:, :, p_i) = BLI_fit.nu_vals;
            end

            Obj.BLI_fit = BLI_fit;
            Obj.nu_pols = BLI(g_nu, 1, Obj.bounds{2}, Obj.order(1, 2), varargin{end});
            %                                       This is i_type_or_x_vals ^^
        end

        function [val, der] = computeWithDer(Obj, x_in)
            [i_nuvals, d_nuvals] = Obj.nu_pols.computeWithDer(x_in(2));

            der = zeros([Obj.op_dims, Obj.in_dims]);
            [val, der(Obj.colons{:}, 1)] = Obj.BLI_fit.computeWithDer(x_in(1), i_nuvals);
            der(Obj.colons{:}, 2) = Obj.BLI_fit.computeWithDer(x_in(1), d_nuvals);
        end

        function der = secondDerivative(Obj, x_in)
            der = 0; % TODO
        end

        function der = firstDerivativeAtPt(Obj, pc_index)
            der = 0; % TODO
        end

        function der = secondDerivativeAtPt(Obj, pc_index)
            der = 0; % TODO
        end

        function save(Obj, filename)
            % Call save on the superclass
            save@Interpolant(Obj, filename);
            BLI_fit = Obj.BLI_fit;
            nu_pols = Obj.nu_pols;

            % Load the remaining class members
            % DEBUG
            % fprintf(2, 'Saving BLI2 Add-ons...\n');
            % END DEBUG
            BLI_fit.save([filename, '.blift' ]); 
            nu_pols.save([filename, '.nupols']);
        end

        function load(Obj, filename)
            % Call load on the superclass
            load@Interpolant(Obj, filename);

            % Load the remaining class members
            % DEBUG
            % fprintf(2, 'Loading BLI2 Add-ons...\n');
            % END DEBUG
            Obj.BLI_fit = BLI();
            Obj.nu_pols = BLI();

            Obj.BLI_fit.load([filename, '.blift' ]);
            Obj.nu_pols.load([filename, '.nupols']);
        end
    end
end
