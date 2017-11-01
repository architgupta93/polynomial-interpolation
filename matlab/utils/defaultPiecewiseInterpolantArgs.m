function args = defaultInterpolantArguments(n_dims)
% function args = defaultInterpolantArguments(n_dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolant Arguments
% Default values to be passed for class instantiation
% Requires 1 input, the number of dimensions for the input
%
% Example:
%   f     = @(x) sin(2*pi*x);
%   args  = defaultInterpolantArguments(1);
%   intrp = BLI(f, args{:});
%
% If no arguments are supplied, arguments are generated for an interpolant with
% a scalar input
%
% Author: Archit Gupta
% Date: October 24, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin == 0)
        n_dims = 1;
    end

    n_pcs   = 5;

    args    = cell(4, 1);
    args{1} = n_dims; % Number of input dimensions. Input can only be a VECTOR

    % The bounds for each dimension
    args{2} = cell(1, n_dims);
    [args{2}{:}] = deal(linspace(-1, 1, n_pcs)');

    % Order of the interpolant
    args{3} = n_pcs * ones(1, n_dims);

    % Category of sample points
    args{4} = 'uniform';
end
