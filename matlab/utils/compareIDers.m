function compareIDers(n_tpts_seed, bounds, df_obj, ip_obj, ip_label, varargin)
% function compareIDers(n_tpts_seed, bounds, df_obj, ip_obj, ip_label, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare the error in derivative both numerically and by plots between a
% function handle and the corresponding interpolant object(s).
% INPUT(s):
%   n_tpts_seed: SCALAR. A seed value for the number of test points PER
%       DIMENSION (the exact number will not be used but something close to it
%       will, this serves as a basic check for code correctness)
%
%   bounds: CELL ARRAY. This should have the same number of cells as the number
%       of inputs that the baseline function handle takes. These values serve
%       as the bounds for test-query points.
%
%   df_obj: FUNCTION HANDLE. Function handle for the derivative of the function
%       to which the polynomial interpolants were fit
%
%   ip_obj: INTERPOLANT. Atleast one interpolant should be supplied for
%       comparing the derivative of the interpolant with the original function.
%
%   ip_label: STRING (optional if only one interpolant is supplied). A label
%       string with which the interpolant will be tagged,
%
%   VARARGIN: A list of other interpolants to compare. This list should be in the form
%       {<ip_obj1>, <ip_label1>, ..., <ip_objn>, <ip_labeln>}
%
% Example:
%   fx        = @(x) sin(2*pi*x);
%   df_dx     = @(x) 2*pi*cos(2*pi*x);
%   bli_ip    = BLI(fx);
%   spline_ip = Spline(fx);
%   lg_ip     = Lagrange(fx);
%   n_tpts    = 1000;
%   bounds    = [-1; 1];
%   compareIDers(n_tpts, bounds, df_dx, bli_ip, 'BLI', spline_ip, 'Spline', ...
%       lg_ip, 'Lagrange');
%
% Author: Archit Gupta, March 09, 2017
% Updated: October 24, 2017 (Moving the plotting to a separate function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    multiple_interpolants = (nargin > 5); 

    if ( nargin < 5 )
        % Just comparing function against one interpolant and label not specified
        ip_label = 'Interpolant';
    end

    % Check that an EVEN number of entries are present in VARARGIN
    if rem(length(varargin), 2)
        error('Expecting an even number of Interpolant/Label entries!')
    end

    n_interpolants = 1 + length(varargin)/2;
    ip_handles     = cell(n_interpolants, 1);
    ip_labels      = cell(n_interpolants, 1);

    % Put the first one in
    ip_handles{1}  = getDervHandle(ip_obj);
    ip_labels{1}   = ip_label;
    df_handle      = getDervHandle(df_obj);

    for ip = 2:n_interpolants
        ip_handles{ip} = getDervHandle(varargin{ip*2-3});
        ip_labels{ip}  = varargin{ip*2 - 2};
    end

    n_in_dims = ip_obj.in_dims;
    n_op_dims = size(df_handle(zeros(n_in_dims, 1)), 1);

    % Generating a set of test points which does not have the same number of
    % points in every dimension. This in itself could be used to find out some
    % basic bugs in the code
    n_test_pts = randi(floor([n_tpts_seed ( n_tpts_seed + 11)]/n_in_dims), ...
        1, n_in_dims);

    g_bounds = zeros(2, n_in_dims);
    for d_i = 1:n_in_dims
        g_bounds(1, d_i) = bounds{d_i}(1);
        g_bounds(2, d_i) = bounds{d_i}(end);    % For piecewise interpolants
    end

    test_obj     = TestingPoints(n_in_dims, n_test_pts, scaleData(g_bounds, 0.99));
    n_test_cases = test_obj.getNPts();

    ip_evals     = cell(n_interpolants, 1);
    ip_evals{:}  = deal(zeros([n_op_dims n_in_dims n_test_pts]));
    df_evals     = zeros([n_op_dims n_in_dims n_test_pts]);

    colons       = cell( size(n_op_dims) );
    [colons{:}]  = deal(':');

    fprintf(2, 'Testing derivatives with %d sample points.\n', n_test_cases); 
    sub = zeros(n_in_dims, 1);
    for t_i = 1:n_test_cases
        [x_in, sub] = test_obj.getPtAt(t_i);
        df_evals(colons{:}, :, sub{:}) = df_handle( x_in );
        for ip = 1:n_interpolants
            ip_evals{ip}(colons{:}, :, sub{:}) = ip_handles{ip}( x_in );
        end
    end

    % Error calculation_
    error_args = ErrorArgs();
    error_args.setDefaults();

    ip_err = zeros(n_interpolants, 1);

    for ip = 1:n_interpolants
        ip_err(ip) = find_error(ip_evals{ip}-df_evals, df_evals, error_args.eps(), error_args.errorMode());
        fprintf(2, '%s Error: %.10f\n', ip_labels{ip}, ip_err(ip));
    end

    s_pts = cell(n_in_dims, 1);
    for d_i = 1:n_in_dims
        s_pts{d_i} = test_obj.getPts(d_i);
    end

    index_to_plot = cell( length(n_op_dims), 1 );
    for p_i = 1 : length( n_op_dims)
        %index_to_plot{p_i} = randi([1 n_op_dims(p_i)]);
        index_to_plot{p_i} = 1;
    end

    if ( n_in_dims == 1)
        fig_handle = plotComparison({}, s_pts, df_evals, ip_evals{:}, 'Baseline', ip_labels{:});
    elseif ( n_in_dims == 2 )
        fprintf(2, 'Plotting dim: %d, df_dx\n', cell2mat(index_to_plot));
        ip_d_dim1 = cell(n_interpolants, 1);
        fig_handle_dim1 = plotComparison({}, s_pts, df_evals(index_to_plot{:}, 1, :, :), ...
            ip_d_dim2{:}, 'Baseline', ip_labels);

        fig_handle_dim2 = plotComparison({}, s_pts, df_evals(index_to_plot{:}, 2, :, :), ...
            ip_d_dim2{:}, 'Baseline', ip_labels);
    else
        fprintf(2, 'Too many dimensions for plotting\n');
    end
end
