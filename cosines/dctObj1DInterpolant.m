classdef dctObj1DInterpolant < Interpolant
% This handle class is solely resposible for handling the "INTERPOLATION" of
% dct-based representation of a function with 1 input argument. It acts as a
% base class for the dctObj1D, which handles extrpolation separately
    properties (Constant = true, Access = public)
        non_zero_threshold = 1e-8;
    end

    properties(Access = private)
        range          = [];
        non_zero_freqs = [];    % Since the signal is expected to have a sparse
                                % representation in the DCT domain, only a few
                                % of the transform coefficient values should
                                % have non-zero values
    end

    methods (Access = protected)
        function non_zero_freqs = getNonZeroFreqs(Obj)
            non_zero_freqs = Obj.non_zero_freqs; 
        end
    end
    
    methods (Access = public)
        function Obj = dctObj1DInterpolant(varargin)
        % function Obj = dctObj1DInterpolant(f_vals, n_in_dims, bounds, order, i_type_or_x_vals);
        % Class constructor
            % Call the parent class constructor first (uncoditionally)
            Obj = Obj@Interpolant(varargin{:});

            % Since we expect the coefficients to be sparse, they are stored as
            % a cell array instead of a Matrix
            Obj.wts            = {};
            Obj.range          = zeros(Obj.in_dims, 1);
            % init_val_dct       = amateur__DCT2(Obj.f_vals', 'MODE_EXTENDED');
            init_val_dct       = dct(Obj.f_vals');
            Obj.non_zero_freqs = cell(Obj.in_dims, 1);

            % Cycle through each input dimension
            for dim_index = 1:Obj.in_dims
                non_zero_indices = find(abs(init_val_dct(:,dim_index)) > ...
                    dctObj1DInterpolant.non_zero_threshold * max(abs(init_val_dct(:,dim_index))));
                if (non_zero_indices(1) == 1)
                    init_val_dct(dim_index, 1) = sqrt(1/2)*init_val_dct(dim_index, 1);
                end
                Obj.wts{dim_index} = sqrt(2/Obj.getNPts(dim_index))*init_val_dct(non_zero_indices, dim_index);
                % Substract 1 to make them start from 1 (and not 2 afer DC)
                non_zero_indices = non_zero_indices - 1;
                Obj.non_zero_freqs{dim_index} = pi*non_zero_indices/Obj.getNPts(dim_index);
                Obj.range(dim_index) = Obj.bounds{dim_index}(2) - Obj.bounds{dim_index}(1);
            end
        end

        % API function calls
        % Interpolation
        function [val, der] = computeWithDer(Obj, xq)
            % Inverse DCT definition (Discrete Version)
            % Xk = 0.5*x0 + sum_{i}(xi*cos(pi*i*(k+0.5)/n_samples)),
            % where Xi is the value of the signal at the index i and xi are the DCT
            % coefficients. This can be generalized to an arbitrary value of x by
            % putting in an appropriate value for i.

            val = zeros([Obj.in_dims 1]);

            for dim_index = 1:Obj.in_dims
                x_scaled_shifted  =  - 0.5 + ( (Obj.getNPts(dim_index)-1) * ...
                    ( xq(dim_index) - Obj.bounds{dim_index}(2) ) / Obj.range(dim_index) );

                val(dim_index, 1) = ...
                    sum( Obj.wts{dim_index} .* ...
                        cos( ...
                            Obj.non_zero_freqs{dim_index} * x_scaled_shifted...
                           ) ...
                       );
                der(dim_index, 1) = - (Obj.getNPts(dim_index) - 1) * ...
                    sum( (Obj.non_zero_freqs{dim_index} .* Obj.wts{dim_index}) .* ...
                        sin( ...
                            Obj.non_zero_freqs{dim_index} * x_scaled_shifted ...
                           ) ...
                       ) / Obj.range(dim_index);
            end
        end

        % Second Derivation
        function der = secondDerivative(Obj, x_inf)
            % TODO: Not defined yet, although this should be easy
            der = 0;
        end
    end
end
