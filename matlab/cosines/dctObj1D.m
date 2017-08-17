classdef dctObj1D < dctObj1DInterpolant
    properties (Access = private)
        extrap_slope = [];
        extrap_value = [];
        extrap_range = [];
    end
    
    methods (Access = public)
        function Obj = dctObj1D(init_val, range_x, extrap_range)
            % init_val is the initial value of the DCT coefficients and range_x is
            % an array of the form [x_min x_max], where x_min and x_max specify the
            % range over which the original function was sampled

            Obj = Obj@dctObj1DInterpolant(init_val, range_x);

            if (nargin == 0)
                % Only need to initialize the quantities that are specific to
                % this class (Not the base-class properties as well. We are
                % happy with the defaults there)
                Obj.extrap_range = [-Inf Inf];
                Obj.extrap_slope = [NaN NaN];
                Obj.extrap_value = [NaN NaN];
            else
                if (nargin < 3)
                    Obj.extrap_range = range_x;
                else
                    Obj.extrap_range = extrap_range;
                end
                
                % This is pretty confusing. The Obj on the right hand side is
                % an instance of the superclass and the one on the left is an
                % instance of the current class.
                n_dims = Obj.getNDims();
                Obj.extrap_value = Inf(n_dims, 2);
                Obj.extrap_slope = Inf(n_dims, 2);

                [Obj.extrap_value(:,1), Obj.extrap_slope(:,1)] = Obj.evaluate(Obj.extrap_range(1));
                [Obj.extrap_value(:,2), Obj.extrap_slope(:,2)] = Obj.evaluate(Obj.extrap_range(2));
            end
        end

        function [yq, dyq_dx] = evaluate(Obj, xq)

            % Inverse DCT definition (Discrete Version)
            % Xk = 0.5*x0 + sum_{i}(xi*cos(pi*i*(k+0.5)/n_samples)),
            % where Xi is the value of the signal at the index i and xi are the DCT
            % coefficients. This can be generalized to an arbitrary value of x by
            % putting in an appropriate value for i.

            if (xq < Obj.extrap_range(:,1))
                yq = Obj.extrap_value(:,1) + Obj.extrap_slope(:,1)*(xq - Obj.extrap_range(1));
                dyq_dx = Obj.extrap_slope(:,1);
            elseif (xq > Obj.extrap_range(2))
                yq = Obj.extrap_value(:,2) + Obj.extrap_slope(:,2)*(xq - Obj.extrap_range(2));
                dyq_dx = Obj.extrap_slope(:,2);
            else
                [yq, dyq_dx] = Obj.interpolate(xq);
            end
        end
    end
end
