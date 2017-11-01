function [f_handle, df_handle, bkpts] = getTestFHandle(n_in_dims, n_op_dims, f_seed)
% function f_handle = getTestFHandle(n_in_dims, n_op_dims, f_seed)
% Author: Archit Gupta, March 03, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Takes number of input dimensions and the number of output dimensions as a
% function and creates a function handle for a function that produces a vector
% ouptut where the size of the output is given by the parameter n_op_dims
% INPUT(s):
%   n_in_dims: [SCALAR] Number of input dimensions
%   n_op_dims: [SCALAR] Number of output dimensions
%   f_seed: Seed for the function to be chosen for test. Options are:
%       - 'c': Continuous (Has a derivative discontinuity)
%       - 'd': Discontinuous (Has a value discontinuity)
%       - 'g': Gaussian (Default. Chosen if an unrecognized 'f_seed' is supplied)
%       - 'o': Oscillatory/Periodic (Sinusoidal function)
%       - 'p': Product Peak (See help f_product_peak)
%       - 'r': Corner Peak (See help f_corner_peak)
%       - 'smooth': one of the smooth functions mentioned above
%
% OUTPUT(s):
%   f_handle: [FUNCTION HANDLE] the function to be tested. This will be one of
%       the functions described in the f_seed description above.
%   df_handle: [FUNCTION HANDLE] Handle for the derivative of the function
%       supplied as f_handle
%   bkpts: [ARRAY] If the supplied function is discontinuous, then bkpts
%       returns the breakpoints, or points of discontinuity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE:
%   1. Checking out BLI interpolation using test function(s)
%   [fun, dfun] = getTestFHandle(1, 1, 'smooth');
%   args        = defaultInterpolantArgs(1);
%   bli_ip      = BLI(fun, args{:});
%   compareIObjs(10, args{2}, fun, bli_ip);
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
        else
            % TODO: Check that f_seed is actually in the available seeds.
        end 
    end

    % Different function handles to try, look up -
    % utils/test_functions in agStuff
    % Choose a function class to test

    c_vec = 4*randn(n_op_dims, n_in_dims);
    w_vec = 1*randn(n_op_dims, n_in_dims);

    switch (f_seed)
        case 'c'
            f_instance = f_continuous(n_in_dims, n_op_dims, [w_vec c_vec]);
        case 'd'
            f_instance = f_discontinuous(n_in_dims, n_op_dims, [w_vec(:, 1:min(2, n_in_dims)) c_vec]);
        case 'g'
            f_instance = f_gaussian(n_in_dims, n_op_dims, [w_vec c_vec]);
        case 'o'
            f_instance = f_oscillatory(n_in_dims, n_op_dims, [w_vec(:, 1) c_vec]);
        case 'r'
            f_instance = f_product_peak(n_in_dims, n_op_dims, [w_vec c_vec]);
        case 'p'
            f_instance = f_corner_peak(n_in_dims, n_op_dims, c_vec);
        otherwise
            f_instance = f_gaussian(n_in_dims, n_op_dims, [w_vec c_vec]);
    end

    % TODO: Print the information about the selected function (or atleast return it)
    bkpts = f_instance.getBreakpoints();
    f_handle = @(x_in) f_instance.evaluate(x_in);
    df_handle = @(x_in) f_instance.diff(x_in);
end
