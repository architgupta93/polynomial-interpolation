classdef ErrorArgs < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% class definition for error calculation arguments for various polynomial
% interpolation experiments.
% Important class members:
%   threshold: Values below threshold * MAX_VALUE are ignore for error
%       calculation
%   err_mode: Let's user select between
%       (MEAN, MAX, MEDIAN) x (REL, ABS) error modes for calculation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = private)
        threshold = [];
        err_mode = [];
    end

    methods (Access = public)
        function args = ErrorArgs(err, th)
            % The user has the option to pass in just one argument 'defaults'
            % which will set the error measurement arguments to their default
            % values
            if (nargin > 1)
                args.threshold = th;
                args.err_mode = err;
            else
                args.threshold = [];
                args.err_mode = [];
                if (nargin > 0)
                    if (strcmp(err, 'defaults'))
                        args.setDefaults();
                    end
                end
            end
        end
        
        function setDefaults(obj)
            obj.threshold = 1e-4;
            obj.err_mode = 'MEAN_{REL}';
        end

        function setMode(obj, mode)
            obj.err_mode = mode;
        end

        function setEps(obj, th)
            obj.threshold = th; 
        end

        function mod = errorMode(obj)
            mod = obj.err_mode;
        end

        function th = eps(obj)
            th = obj.threshold;
        end
    end
end
