classdef f_discontinuous < HDFunction
    properties (Access = protected)
        c_vec = [];
        w_vec = [];
        n_elems_to_chk = 0;
    end

    methods (Access = public)
        function Obj = f_discontinuous(n_in, n_out, parms)
            Obj = Obj@HDFunction(n_in, n_out, parms);
            Obj.n_elems_to_chk = min(2, Obj.n_in);
            if ( size(Obj.parms, 2) ~= Obj.n_elems_to_chk+Obj.n_in )
                fprintf('ERROR: Mismatch in the parameter/input dimensions\n');
                fprintf('Expecting [w1 w2 c_vec]\n');
                return;
            elseif ( size(Obj.parms, 1) ~= Obj.n_out )
                fprintf('ERROR: Mismatch in the parameter/output dimensions\n');
                return;
            end
            Obj.w_vec = Obj.parms(:, 1:Obj.n_elems_to_chk);
            Obj.c_vec = Obj.parms(:, Obj.n_elems_to_chk+1:end);
        end

        function y_out = evaluate(Obj, x_in)
            Obj.checkDims(x_in);
            y_out = exp(Obj.c_vec * x_in) .* prod(Obj.w_vec(:, 1:Obj.n_elems_to_chk) < x_in(1:Obj.n_elems_to_chk)', 2);
        end

        function der = diff(Obj, x_in)
            Obj.checkDims(x_in);
            der = times( Obj.c_vec, exp(Obj.c_vec * x_in) .* prod(Obj.w_vec(:, 1:Obj.n_elems_to_chk) < ...
                x_in(1:Obj.n_elems_to_chk)', 2) );
        end
    end
end
