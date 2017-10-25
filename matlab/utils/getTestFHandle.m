function [f_handle, df_handle, bkpts] = getTestFHandle(n_in_dims, n_op_dims, f_seed)
% function f_handle = getTestFHandle(n_in_dims, n_op_dims, f_seed)
% Author: Archit Gupta, March 03, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Takes number of input dimensions and the number of output dimensions as a function and creates a function handle for a
% function that produces a vector ouptut where the size of the output is given by the parameter n_op_dims
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    narginchk(2, 3);
    N_FUNS = 6;
    SEEDS = ['c', 'd', 'g', 'o', 'r', 'p'];
    if (nargin < 3)
        i_seed = randi([1, N_FUNS], 1, 1);
        f_seed = SEEDS(i_seed);
    else
        if ( strcmp(f_seed, 'smooth') ) % Limit the choice of functions to smooth functions
            % These functions include 'g', 'o', 'r' ,'p'
            i_seed = randi([3, N_FUNS], 1, 1);
            f_seed = SEEDS(i_seed);
        end 
    end

    % Different function handles to try, look up -
    % utils/test_functions in agStuff
    % Choose a function class to test

    c_vec = 4*randn(n_op_dims, n_in_dims);
    w_vec = 1*randn(n_op_dims, n_in_dims);

    switch (f_seed)
        case 'c'
            f_instance = f_continuous(n_in_dims, n_op_dims, [w_vec c_vec])
        case 'd'
            f_instance = f_discontinuous(n_in_dims, n_op_dims, [w_vec(:, 1:min(2, n_in_dims)) c_vec])
        case 'g'
            f_instance = f_gaussian(n_in_dims, n_op_dims, [w_vec c_vec])
        case 'o'
            f_instance = f_oscillatory(n_in_dims, n_op_dims, [w_vec(:, 1) c_vec])
        case 'r'
            f_instance = f_product_peak(n_in_dims, n_op_dims, [w_vec c_vec])
        case 'p'
            f_instance = f_corner_peak(n_in_dims, n_op_dims, c_vec)
        otherwise
            f_instance = f_gaussian(n_in_dims, n_op_dims, [w_vec c_vec])
    end

    bkpts = f_instance.getBreakpoints();
    f_handle = @(x_in) f_instance.evaluate(x_in);
    df_handle = @(x_in) f_instance.diff(x_in);
end
