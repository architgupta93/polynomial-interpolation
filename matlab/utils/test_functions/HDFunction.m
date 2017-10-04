classdef HDFunction < handle
    properties (Access = protected)
        n_in = 0;
        n_out = 0;
        parms = [];
    end

    methods (Access = protected)
        function checkDims(Obj, x_in)
            assert( size(x_in, 1) == Obj.n_in, 'ERROR: x_in incorrectly sized');
        end
    end
    
    methods (Access = public)
        function Obj = HDFunction(n_in, n_out, parms)
            Obj.n_in = n_in;
            Obj.n_out = n_out; 
            Obj.parms = parms;
        end
    end

    methods (Access =  public) % Abstract function definitions from MATLAB
        % Every implementation of this class must implement the () operator, or
        % in case of Matlab, the subsindex function
        function y_out = evaluate(Obj, x_in)
            y_out = [];
        end

        function der = diff(Obj, x_in)
            der = [];
        end

        function bkpts = getBreakpoints(Obj)
            bkpts = [];
        end
    end
end
