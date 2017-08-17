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

    if ( dim > length(size(v_in)) )
        error('Supplied DIM exceeds data dimensions');
    end

    n_elems_in_dim = size(v_in, dim);
    d_mean = mean(v_in, dim);

    dim_array = ones([1 length(size(v_in))]);
    dim_array(dim) = n_elems_in_dim;
    c_data = v_in - repmat(d_mean, dim_array);
    v_out = (scale * c_data) + repmat(d_mean, dim_array);
end
