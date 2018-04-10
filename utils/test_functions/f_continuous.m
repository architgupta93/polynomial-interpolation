classdef f_continuous < HDFunction
    properties (Access = protected)
        c_vec = [];
        w_vec = [];
    end

    methods (Access = public)
        function Obj = f_continuous(n_in, n_out, parms)
            Obj = Obj@HDFunction(n_in, n_out, parms);
            if ( size(Obj.parms, 2) ~= 2*Obj.n_in )
                fprintf('ERROR: Mismatch in the parameter/input dimensions\n');
                fprintf('Expecting [w_vec c_vec]\n');
                return;
            elseif ( size(Obj.parms, 1) ~= Obj.n_out )
                fprintf('ERROR: Mismatch in the parameter/output dimensions\n');
                return;
            end
            Obj.w_vec = Obj.parms(:, 1:Obj.n_in);
            Obj.c_vec = Obj.parms(:, 1+Obj.n_in:end);
        end

        function y_out = evaluate(Obj, x_in)
            Obj.checkDims(x_in);
            y_out = exp(-sum(times(Obj.c_vec, abs(minus(x_in', Obj.w_vec))), 2));
        end

        function der = diff(Obj, x_in)
            Obj.checkDims(x_in);
            der = times( -times(Obj.c_vec, sign(minus(x_in', Obj.w_vec))), ...
                exp(-sum(times(Obj.c_vec, abs(minus(x_in', Obj.w_vec))), 2)) );
        end

        function bkpts = getBreakpoints(Obj)
            bkpts = Obj.w_vec;
        end
    end
end
