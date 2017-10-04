classdef f_corner_peak < HDFunction
    properties (Access = protected)
        c_vec = [];
    end

    methods (Access = public)
        function Obj = f_corner_peak(n_in, n_out, parms)
            Obj = Obj@HDFunction(n_in, n_out, parms);
            if ( size(Obj.parms, 2) ~= Obj.n_in )
                fprintf('ERROR: Mismatch in the parameter/input dimensions\n');
                fprintf('Expecting [c_vec]\n');
                return;
            elseif ( size(Obj.parms, 1) ~= Obj.n_out )
                fprintf('ERROR: Mismatch in the parameter/output dimensions\n');
                return;
            end
            Obj.c_vec = Obj.parms;    
        end

        function y_out = evaluate(Obj, x_in)
            Obj.checkDims(x_in);
            y_out = ( 10 + Obj.c_vec * x_in ).^(-(1+Obj.n_in));
        end

        function der = diff(Obj, x_in)
            Obj.checkDims(x_in);
            der = (-(1 + Obj.n_in)) * times( Obj.c_vec, ( 10 + Obj.c_vec * x_in ).^(-(2+Obj.n_in)) );
        end
    end
end
