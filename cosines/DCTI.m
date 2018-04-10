classdef DCTI < dctObj1DInterpolant
    properties (Access = private)
        extrap_value = [];

        % Since DCT ends up being very sloppy towards the boundary, the
        % extrapolation range could also overlap with the range in which we
        % sampled the function values
        extrap_range = [];
    end
    
    methods (Access = public)
        function Obj = DCTI(varargin)
        % function Obj = dctObj1D(f_vals, n_in_dims, bounds, order, i_type_or_x_vals);
        % Class constructor
            % Call the parent class constructor
            Obj = Obj@dctObj1DInterpolant(varargin{:});

            if (isempty(varargin))
                % Only need to initialize the quantities that are specific to
                % this class (Not the base-class properties as well. We are
                % happy with the defaults there)
                Obj.extrap_range = [-Inf Inf];
                Obj.extrap_value = [NaN NaN];
            else
                if (nargin > 5)
                    % If extrapolation range is supplied externally, use that
                    Obj.extrap_range = varargin{6}{1};
                else
                    % Otherwise, use the boundary value
                    % HACK, Fix later for  multiple inputs... TODO
                    Obj.extrap_range = Obj.bounds{1};
                end
                
                % This is pretty confusing. The Obj on the right hand side is
                % an instance of the superclass and the one on the left is an
                % instance of the current class.
                Obj.extrap_value = Inf(Obj.in_dims, 2);
                Obj.extrap_slope = Inf(Obj.in_dims, 2);

                [Obj.extrap_value(:,1), Obj.extrap_slope(:,1)] = Obj.computeWithDer(Obj.extrap_range(1));
                [Obj.extrap_value(:,2), Obj.extrap_slope(:,2)] = Obj.computeWithDer(Obj.extrap_range(2));
            end
        end

        % Adding the implementation of API function (extrapolation included)
        function [val, der] = computeWithDer(Obj, xq)

            % Inverse DCT definition (Discrete Version)
            % Xk = 0.5*x0 + sum_{i}(xi*cos(pi*i*(k+0.5)/n_samples)),
            % where Xi is the value of the signal at the index i and xi are the DCT
            % coefficients. This can be generalized to an arbitrary value of x by
            % putting in an appropriate value for i.

            if (xq < Obj.extrap_range(:,1))
                val = Obj.extrap_value(:,1) + Obj.extrap_slope(:,1)*(xq - Obj.extrap_range(1));
                der = Obj.extrap_slope(:,1);
            elseif (xq > Obj.extrap_range(2))
                val = Obj.extrap_value(:,2) + Obj.extrap_slope(:,2)*(xq - Obj.extrap_range(2));
                der = Obj.extrap_slope(:,2);
            else
                [val, der] = computeWithDer@dctObj1DInterpolant(Obj, xq);
            end
        end

        % Second derivative will also be different with extrapolation
        function der = secondDerivative(Obj, x_in)
            val = [];
            der = [];
        end
    end
end
