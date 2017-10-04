classdef f_oscillatory < HDFunction
    properties (Access = protected)
        w = 0;
        c_vec = [];
    end

    methods (Access = public)
        function Obj = f_oscillatory(n_in, n_out, parms)
            Obj = Obj@HDFunction(n_in, n_out, parms);
            if ( size(Obj.parms, 2) ~= 1+Obj.n_in )
                fprintf('ERROR: Mismatch in the parameter/input dimensions\n');
                fprintf('Expecting [w c_vec]\n');
                return;
            elseif ( size(Obj.parms, 1) ~= Obj.n_out )
                fprintf('ERROR: Mismatch in the parameter/output dimensions\n');
                return;
            end

            Obj.w = Obj.parms(:, 1);
            Obj.c_vec = Obj.parms(:, 2:end);  % Store as a row vector and avoid
                                            % transpose everytime required
        end

        function y_out = evaluate(Obj, x_in)
            Obj.checkDims(x_in);
            y_out = cos( 2*pi*Obj.w + Obj.c_vec * x_in );
        end

        function der = diff(Obj, x_in)
            Obj.checkDims(x_in);
            der = times( Obj.c_vec, -sin( 2*pi*Obj.w + Obj.c_vec * x_in ) );
        end
    end
end
