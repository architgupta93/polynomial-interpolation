classdef dctObj2DInterpolant < handle
    properties (Constant, Access = public)
        non_zero_threshold = 1e-8;
    end

    properties(Access = private)
        % In general, one should be using tensors for the following quantities
        % but let's see how far we can go with this simple hacked version
        x_min = [];
        x_max = [];
        y_min = [];
        y_max = [];
    end

    properties(SetAccess = immutable)
        n_dims = [];
        n_samples = [];
        range_x = [];
        range_y = [];
        coeffs = [];
        non_zero_freqs_x = [];
        non_zero_freqs_y = [];
    end

    methods(Access = protected)
        function n_dims = getNDims(Obj)
            n_dims = Obj.n_dims;
        end

        function n_samples = getNSamples(Obj)
            n_samples = Obj.n_samples;
        end

        function [freqs_x freqs_y] = getNonZeroFreqs(Obj)
            freqs_x = Obj.non_zero_freqs_x;
            freqs_y = Obj.non_zero_freqs_y;
        end

        function coeffs = getCoeffs(Obj)
            coeffs = Obj.coeffs;
        end

        function x_index = findXIndex(Obj, x_val)
            x_index = (Obj.n_samples(1)-1).*(x_val - Obj.x_min)./Obj.range_x;
        end

        function y_index = findYIndex(Obj, y_val)
            y_index = (Obj.n_samples(2)-1).*(y_val - Obj.y_min)./Obj.range_y;
        end
    end

    methods(Access = public)
        function Obj = dctObj2DInterpolant(init_val, range_x, range_y)
            if (nargin == 0)
                % Set a few default values for safety
                Obj.n_samples = [0 0];
                Obj.n_dims = 0;
                Obj.x_min = -Inf;
                Obj.x_max = Inf;
                Obj.y_min = -Inf;
                Obj.y_max = Inf;
                Obj.coeffs = 0;
                Obj.non_zero_freqs_x = nan;
                Obj.non_zero_freqs_y = nan;
            else
                Obj.n_samples(1) = size(init_val, 1);
                Obj.n_samples(2) = size(init_val, 2);
                Obj.n_dims = size(init_val, 3);

                % TODO: The next few lines will probably come back and bite us
                % if we don't get them right. Remember that for the MOSFET data,
                % the order of magnitude of the fi/fe functions is very
                % different from qe/qi functions. So while ignoring DCT values,
                % we need to make sure that the comparison is only made between
                % all the values corresponding to the same dimension

                % It is important to maintain the 2D structre for reconstruction
                Obj.non_zero_freqs_x = cell(Obj.n_dims, 1);
                Obj.non_zero_freqs_y = cell(Obj.n_dims, 1);

                Obj.coeffs = cell(Obj.n_dims, 1);
                for dim_index = 1:Obj.n_dims
                    dct_coeffs_in_this_dim = dct2(init_val(:,:,dim_index));
                    max_val_in_this_dim = max(max(abs(dct_coeffs_in_this_dim)));
                    [non_zero_indices_x, non_zero_indices_y] = ...
                        find( ...
                                abs(dct_coeffs_in_this_dim) > ...
                                dctObj2DInterpolant.non_zero_threshold * max_val_in_this_dim ...
                            );
                    % Rescaling the coefficients (Quite painful as the DC Values
                    % and the others have to be handle separately)
                    coefficient_scaling_factor = (2/sqrt(prod(Obj.n_samples)))*ones(size(non_zero_indices_x));  % Same as size(non_zero_indices_y)
                    dc_vals_in_x_dim = find(non_zero_indices_x == 1);
                    dc_vals_in_y_dim = find(non_zero_indices_y == 1);
                    
                    if (length(dc_vals_in_x_dim))
                        coefficient_scaling_factor(dc_vals_in_x_dim) = sqrt(1/2)*coefficient_scaling_factor(dc_vals_in_x_dim);
                    end
                    if (length(dc_vals_in_y_dim))
                        coefficient_scaling_factor(dc_vals_in_y_dim) = sqrt(1/2)*coefficient_scaling_factor(dc_vals_in_y_dim);
                    end

                    % Handling the coefficients on a dimension-by-dimension basis
                    Obj.coeffs{dim_index} = ...
                        coefficient_scaling_factor .* ...
                        dct_coeffs_in_this_dim( ...
                                non_zero_indices_x + ...
                                Obj.n_samples(1)*(non_zero_indices_y  - 1) ...
                            );

                    Obj.non_zero_freqs_x{dim_index} = pi*(non_zero_indices_x - 1) / Obj.n_samples(1);
                    Obj.non_zero_freqs_y{dim_index} = pi*(non_zero_indices_y - 1) / Obj.n_samples(2);
                end

                if ( (range_x(1) == range_x(2)) | (range_y(1) == range_y(2)) )
                    fprintf('ERROR: Invalid data range supplied!\n');
                end

                Obj.x_min = range_x(1);
                Obj.x_max = range_x(2);
                Obj.y_min = range_y(1);
                Obj.y_max = range_y(2);
                Obj.range_x = Obj.x_max - Obj.x_min;
                Obj.range_y = Obj.y_max - Obj.y_min;
                % Calculating compression gained
                original_dimensions = prod(size(init_val));
                n_stored_coeffs = 0;
                for dim_index = 1:Obj.n_dims
                    n_stored_coeffs = n_stored_coeffs + prod(size(Obj.coeffs{dim_index}));
                end
                fprintf('Using DCT interpolant: Obtained a compression of %.1fX\n', original_dimensions/n_stored_coeffs);
            end
        end

        function [vq_calc, dvq_dx_calc, dvq_dy_calc] = interpolate(Obj, xq, yq)
            % Even in 2 dimensions, DCT is simply an extension of DCT (type 2),
            % applied row-wise (and then column-wise on the resulting
            % coefficients). It is very similar to the idea behind tensor
            % product splines. Just implementing DCT (Type 3) in the same way
            % should do the job for inverting it

            vq_calc = zeros([Obj.n_dims 1]);
            dvq_dx_calc = zeros([Obj.n_dims 1]);
            dvq_dy_calc = zeros([Obj.n_dims 1]);

            index_x = Obj.findXIndex(xq);
            index_y = Obj.findYIndex(yq);
            for dim_index = 1:Obj.n_dims
                vq_calc(dim_index, 1) = ...
                    sum ( ...
                    sum ( Obj.coeffs{dim_index} .* ...
                            cos ( ...
                                    Obj.non_zero_freqs_x{dim_index} * (index_x + 0.5) ...
                                ) .* ...
                            cos ( ...
                                    Obj.non_zero_freqs_y{dim_index} * (index_y + 0.5) ...
                                ) ... 
                        )); 

                if ( nargout > 1)   % Its too much pain to be taken unnecessarily
                    dvq_dy_calc(dim_index, 1) = - (Obj.n_samples(2) - 1) * ...
                        sum ( ...
                        sum ( Obj.non_zero_freqs_y{dim_index} .* Obj.coeffs{dim_index} .* ...
                                cos ( ...
                                        Obj.non_zero_freqs_x{dim_index} * (index_x + 0.5) ...
                                    ) .* ...
                                sin ( ...
                                        Obj.non_zero_freqs_y{dim_index} * (index_y + 0.5) ...
                                    )  ... 
                            )) / Obj.range_y; 

                    dvq_dx_calc(dim_index, 1) = - (Obj.n_samples(1) - 1) * ...
                        sum ( ...
                        sum ( Obj.non_zero_freqs_x{dim_index} .* Obj.coeffs{dim_index} .* ...
                                sin ( ...
                                        Obj.non_zero_freqs_x{dim_index} * (index_x + 0.5) ...
                                    ) .* ...
                                cos ( ...
                                        Obj.non_zero_freqs_y{dim_index} * (index_y + 0.5) ...
                                    ) ...  
                            )) / Obj.range_x; 
                end
            end
        end
    end
end
