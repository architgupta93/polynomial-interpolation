classdef UniformPoints < InterpolationPoints
    properties (SetAccess = protected)
        step_size = [];
    end

    methods (Access = protected)
        function pts = generate(Obj, order)
            pts = cell(Obj.n_in_dims, 1);
            n_pts = 2.^order + 1;
            for d_i = 1:Obj.n_in_dims
                pts{d_i} = linspace( 1, -1, n_pts(d_i) )';
            end
        end
    end

    methods (Access = public)
        % Class constructor
        function Obj = UniformPoints(n_in_dims, order, bounds)
            Obj = Obj@InterpolationPoints(n_in_dims, order, bounds);
            pts = Obj.generate(Obj.order);
            Obj.recenterAndSave(pts);

            % Some class properties required for generating BLI Weights
            Obj.step_size = zeros(1, n_in_dims);
            for d_i = 1:n_in_dims
                Obj.step_size(1, d_i) = Obj.getPtAt(2, d_i) - Obj.getPtAt(1, d_i);
            end
        end

        % Get weights for BLI
        function wts = getBLIWeights(obj)
            % BLI weights (for index i) are obtained by the following formula
            %
            %               wj = 1 / prod_(k ~= j) (xj - xk)
            %
            % If the points are uniformly distributed (as in this case), This
            % can be reduced to (just simple counting, assuming that the indices
            % for j vary in [1, 2, ..., n]):
            %               
            %               wj = 1 / ( (j-1)!(n-j)! * (dx ^ (n-1)) )
            %
            % Note that this cannot be done explicitly (using factorial formulae
            % for (j-1) and (n-j) as Matlab's factorial(200) is Inf and it
            % generally seems like a bad idea

            wts = cell(obj.n_in_dims, 1);
            for d_i = 1:obj.n_in_dims
                s_factorials = ones(obj.n_pts(d_i), 1);
                s_factorials(2,1) = obj.step_size(d_i);
                for f_i = 2:obj.n_pts(d_i)-1
                    s_factorials(f_i+1,1) = s_factorials(f_i,1) * f_i * obj.step_size(d_i);
                end

                wts{d_i} = zeros(obj.n_pts(d_i), 1);
                one = 1;
                for f_i = 1:obj.n_pts(d_i)
                    one = -1 * one;
                    wts{d_i}(f_i) = one / ( s_factorials(f_i) * ...
                        s_factorials(1+obj.n_pts(d_i)-f_i) );
                end
            end
        end
    end
end
