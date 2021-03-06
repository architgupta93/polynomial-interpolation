classdef BLI2 < Interpolant
    properties (Access = protected)
        BLI_fit   = BLI();
        coeff_fit = BLI();

        bli_fit_ext   = '.blift';
        coeff_fit_ext = '.coeff';
    end

    methods (Access = public)
        function Obj = BLI2(varargin)
        % function Obj = BLI2(f_vals, in_dims, bounds, order, i_type_or_x_vals) 
            Obj = Obj@Interpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            % Express the function values as a BLI expression in the first dimensioN
            coeffs_1d = zeros([Obj.op_dims, Obj.i_pts.n_pts(1), Obj.i_pts.n_pts(2)]);
            
            for p_i = 1:Obj.i_pts.n_pts(2)
                BLI_fit = BLI(Obj.f_vals(Obj.colons{:}, :, p_i), 1, {Obj.bounds{1}}, Obj.order(1, 1), varargin{end});
                %                                                           This is i_type_or_x_vals ^^
                coeffs_1d(Obj.colons{:}, :, p_i) = BLI_fit.coeffs;
            end

            Obj.BLI_fit   = BLI_fit;
            Obj.coeff_fit = BLI(coeffs_1d, 1, {Obj.bounds{2}}, Obj.order(1, 2), varargin{end});
            %                                       This is i_type_or_x_vals ^^
        end

        function [val, der] = computeWithDer(Obj, x_in)
            if (nargout > 1)
                der                          = zeros([Obj.op_dims, Obj.in_dims]);
                [i_nuvals, d_nuvals]         = Obj.coeff_fit.computeWithDer(x_in(2));
                [val, der(Obj.colons{:}, 1)] = Obj.BLI_fit.computeWithDer(x_in(1), i_nuvals);
                der(Obj.colons{:}, 2)        = Obj.BLI_fit.computeWithDer(x_in(1), d_nuvals);
                return;
            end

            val = Obj.BLI_fit.computeWithDer(x_in(1), ...
                Obj.coeff_fit.computeWithDer(x_in(2)));
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
            BLI_fit   = Obj.BLI_fit;
            coeff_fit = Obj.coeff_fit;

            % Load the remaining class members
            % DEBUG
            % fprintf(2, 'Saving BLI2 Add-ons...\n');
            % END DEBUG
            BLI_fit.save([filename, Obj.bli_fit_ext]); 
            coeff_fit.save([filename, Obj.coeff_fit_ext]);
        end

        function load(Obj, filename)
            % Call load on the superclass
            load@Interpolant(Obj, filename);

            % Load the remaining class members
            % DEBUG
            % fprintf(2, 'Loading BLI2 Add-ons...\n');
            % END DEBUG
            Obj.BLI_fit = BLI();
            Obj.coeff_fit = BLI();

            Obj.BLI_fit.load([filename, Obj.bli_fit_ext]);
            Obj.coeff_fit.load([filename, Obj.coeff_fit_ext]);
        end
    end
end
