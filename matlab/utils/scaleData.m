function v_out = scaleData(v_in, scale, dim)
% function v_out = shrinkData(v_in, scale)
% Author: Archit Gupta (Feb 28, 2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function takes in a vector/matrix/tensor and squeezes the data inside it by a constant factor given by the factor
% "scale". This is done by centering the data first, scaling it and then moving it b ack to the original center. The
% user can optionally supply a "dim" along which the data should be centered for this process
%
% INPUTS:
%   v_in: Data that need to be scaled
%   scale: Factor by which the data needs to be scaled
%   dim: (optional) dimension along which the data should be centered
%
% OUTPUS:
%   v_out: Rescaled data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    narginchk(2,3);

    if (nargin < 3)
        dim = 1;
    end

    % If v_in is a cell array, perform this operation only on the matrix along 'dim'.
    % Since this removes the original functionality of DIM, for column
    % vectors, dim will be set to 1, otherwise, dim=2 will be used whenever
    % v_in is a cell array.
    cell_to_mat_conversion = false;
    if (iscell(v_in))
        v_in_cell = v_in;
        v_in      = v_in_cell{dim};
        d_cell    = dim;

        % Change DIM so that row/column vectors are handled properly
        % Applicable only when v_in is a cell array.
        if iscolumn(v_in)
            dim = 1;
        else
            dim = 2;
        end
        cell_to_mat_conversion = true;
    end

    if ( dim > length(size(v_in)) )
        error('Supplied DIM exceeds data dimensions');
    end

    n_elems_in_dim = size(v_in, dim);
    d_mean = mean(v_in, dim);

    dim_array = ones([1 length(size(v_in))]);
    dim_array(dim) = n_elems_in_dim;
    c_data = v_in - repmat(d_mean, dim_array);
    v_out = (scale * c_data) + repmat(d_mean, dim_array);

    if (cell_to_mat_conversion)
        % Replace the scale dimension in the cell array.
        v_in_cell{d_cell} = v_out;
        v_out = v_in_cell;
    end
end
