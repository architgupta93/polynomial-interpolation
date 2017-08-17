classdef dctObj1DInterpolant < handle
% This handle class is solely resposible for handling the "INTERPOLATION" of
% dct-based representation of a function with 1 input argument. It aacts as a
% base class for the dctObj1D, which handles extrpolation separately
    properties (Constant = true, Access = public)
        non_zero_threshold = 1e-8;
    end

    properties(Access = private)
        n_dims = [];    % Number of dimensions in the function
        n_samples = []; % Length of discrete input signal / cosine-coefficient
                        % vector
        coeffs = [];    % Non-zero DCT Coefficients
        non_zero_freqs = [];    % Since the signal is expected to have a sparse
                                % representation in the DCT domain, only a few
                                % of the transform coefficient values should
                                % have non-zero values
        x_min = [];
        x_max = [];     % x_min and x_max will be used to find the indices of
                        % an arbitrary point. If the point does not lie at a
                        % sampling point, this index will not be an integer
                        % either 
        range = [];
    end

    methods (Access = protected)
        function n_dims = getNDims(Obj)
            n_dims = Obj.n_dims;
        end

        function non_zero_freqs = getNonZeroFreqs(Obj)
            non_zero_freqs = Obj.non_zero_freqs; 
        end

        function n_samples = getNSamples(Obj)
            n_samples = Obj.n_samples; 
        end

        function index = findIndex(Obj, x_val)
            % Works for a vector x_min/x_max/n_samples as well, i.e., for a
            % multi-dimensional function (all dimensions evaluated at a scalar
            % x_val
            index = (Obj.n_samples-1).*(x_val - Obj.x_min)./Obj.range;
        end

        function coeffs = getCoeffs(Obj)
            coeffs = Obj.coeffs;
        end
    end
    
    methods (Access = public)
        function Obj = dctObj1DInterpolant(init_val, range_x)
            % init_val is the initial value of the DCT coefficients and range_x is
            % an array of the form [x_min x_max], where x_min and x_max specify the
            % range over which the original function was sampled

            if (nargin == 0)
                Obj.n_samples = 0;
                Obj.n_dims = 0;
                Obj.x_min = -Inf;
                Obj.x_max = Inf;
                Obj.coeffs = 0;
                Obj.non_zero_freqs = nan;
            else
                Obj.n_samples = size(init_val, 1);
                Obj.n_dims = size(init_val, 2);
                init_val_dct = dct(init_val);
                Obj.non_zero_freqs = cell(Obj.n_dims, 1);
                Obj.coeffs = cell(Obj.n_dims, 1);

                for dim_index = 1:Obj.n_dims
                    non_zero_indices = find(abs(init_val_dct(:,dim_index)) > ...
                        dctObj1DInterpolant.non_zero_threshold * max(abs(init_val_dct(:,dim_index))));
                    if (non_zero_indices(1) == 1)
                        init_val_dct(1,dim_index) = sqrt(1/2)*init_val_dct(1,dim_index);
                    end
                    Obj.coeffs{dim_index} = sqrt(2/Obj.n_samples)*init_val_dct(non_zero_indices, dim_index);
                    % Substract 1 to make them start from 1 (and not 2 afer DC)
                    non_zero_indices = non_zero_indices - 1;
                    Obj.non_zero_freqs{dim_index} = pi*non_zero_indices/Obj.n_samples;
                end

                if (range_x(1) == range_x(2))
                    fprintf('ERROR: Invalid data range supplied!\n');
                end

                Obj.x_min = range_x(1);
                Obj.x_max = range_x(2);
                Obj.range = Obj.x_max - Obj.x_min;
            end
        end

        function [yq_calc, dyq_dx_calc] = interpolate(Obj, xq)

            % Inverse DCT definition (Discrete Version)
            % Xk = 0.5*x0 + sum_{i}(xi*cos(pi*i*(k+0.5)/n_samples)),
            % where Xi is the value of the signal at the index i and xi are the DCT
            % coefficients. This can be generalized to an arbitrary value of x by
            % putting in an appropriate value for i.

            yq_calc = zeros([Obj.n_dims 1]);
            index_for_all_dims = Obj.findIndex(xq);

            for dim_index = 1:Obj.n_dims
                yq_calc(dim_index, 1) = ...
                    sum( Obj.coeffs{dim_index} .* ...
                        cos( ...
                            Obj.non_zero_freqs{dim_index} * (index_for_all_dims + 0.5) ...
                           ) ...
                       );
                dyq_dx_calc(dim_index, 1) = - (Obj.n_samples - 1) * ...
                    sum( (Obj.non_zero_freqs{dim_index} .* Obj.coeffs{dim_index}) .* ...
                        sin( ...
                            Obj.non_zero_freqs{dim_index} * (index_for_all_dims + 0.5) ...
                           ) ...
                       ) / Obj.range;
            end
        end
    end
end
