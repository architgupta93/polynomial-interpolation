function t_prod = tensorProduct(in_array)
% function t_prod = tensor_product(in_array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes a tensor product of several input vectors. All the input
% vectors should be put in a cell-array, called in_array.
%
%   INPUT:
%       in_array: Input cell array that contains all the vectors that need to be
%       multiplied together.
%
%   OUTPUT:
%       t_prod: Tensor product (all possible multiplications)
%
% Author: Archit Gupta (Feb 16, 2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Helper functions below for performing permutations on the dimensions of
    % the various vectors in the tensor product operation above

    n_elems = numel(in_array);

    t_prod = in_array{1};
    for elem = 2:n_elems
        e_t = permute( in_array{elem}, circshift( 1:(g_dims(t_prod) + ...
            g_dims(in_array{elem})), [0, g_dims(in_array{elem})] ) );
        t_prod = bsxfun( @times, t_prod, e_t);
    end
end

