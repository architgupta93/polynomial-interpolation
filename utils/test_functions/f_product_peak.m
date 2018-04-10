classdef f_product_peak < HDFunction
    properties (Access = protected)
        c_vec = [];
        w_vec = [];
    end

    methods (Access = public)
        function Obj = f_product_peak(n_in, n_out, parms)
            Obj = Obj@HDFunction(n_in, n_out, parms);
            if ( size(Obj.parms, 2) ~= 2*Obj.n_in )
                fprintf('ERROR: Mismatch in the parameter/input dimensions\n');
                fprintf('Expecting [w_vec c_vec]\n');
                return;
            elseif ( size(Obj.parms, 1) ~= Obj.n_out )
                fprintf('ERROR: Mismatch in the parameter/output dimensions\n');
                return;
            end
            Obj.w_vec = 2 + Obj.parms(:, 1:Obj.n_in);
            Obj.c_vec = 10 + Obj.parms(:, 1+Obj.n_in:end);
        end

        function y_out = evaluate(Obj, x_in)
            Obj.checkDims(x_in);
            y_out = 1./prod( minus(Obj.c_vec.^2, (minus(x_in', Obj.w_vec)).^2), 2); 
        end

        function der = diff(Obj, x_in)
            Obj.checkDims(x_in);
            der = rdivide( 2 * (minus(x_in', Obj.w_vec)), times(minus(Obj.c_vec.^2, (minus(x_in', Obj.w_vec)).^2), ...
                prod( minus(Obj.c_vec.^2, (minus(x_in', Obj.w_vec)).^2), 2)) );
        end
    end
end
