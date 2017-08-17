classdef ChebyshevPoints < InterpolationPoints
    properties (Access = protected)
        cheb_type = 0;
        roots_or_extremas = false;
        grid_type_dense = true;
                    % The last option (grid_type_dense) implies the use of a
                    % regular, dense, tensor product grid over n-dimensions
    end

    methods (Access = public) % Abstract functions from MATLAB
        % Abstract function for generating new "n_pts" data points 
        function [pts, n_pts] = generate(Obj, order);
            pts = [];
            n_pts = [];
        end
    end

    methods (Access = public)
        function Obj = ChebyshevPoints(n_dims, order, bounds)
            Obj = Obj@InterpolationPoints(n_dims, order, bounds);
            Obj.grid_type_dense = true;
        end
    end
end
