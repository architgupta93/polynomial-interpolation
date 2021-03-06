function [est_err, speedup, err_plt] = compareIObjs(n_tpts_seed, bounds, ...
    f_obj, varargin)
% function [est_err, speedup] = compareIObjs(n_tpts_seed, bounds, ...
%    f_obj, varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare the error both numerically and by plots between a function handle and
% the corresponding interpolant object(s).
% INPUT(s):
%   n_tpts_seed: SCALAR. A seed value for the number of test points PER
%       DIMENSION (the exact number will not be used but something close to it
%       will, this serves as a basic check for code correctness)
%
%   bounds: CELL ARRAY. This should have the same number of cells as the number
%       of inputs that the baseline function handle takes. These values serve
%       as the bounds for test-query points.
%
%   f_obj: FUNCTION HANDLE. Function handle for evaluating the function to
%       which the polynomial interpolants were fit
%
%   VARARGIN: A list of other interpolants to compare. This list should be in the form
%       {<ip_obj1>, <ip_label1>, ..., <ip_objn>, <ip_labeln>}, where ip_obj<i> and
%       ip_label<i> are:
%
%   ip_obj: INTERPOLANT. Atleast one interpolant should be supplied for
%       comparing the derivative of the interpolant with the original function.
%
%   ip_label: STRING (optional if only one interpolant is supplied). A label
%       string with which the interpolant will be tagged,
%
% OUTPUT(s):
%   EST_ERR: estimated error between the function nad the interpolant objects
%   SPEEDUP: Relative time taken by the interpolants while having the
%       function's time as baseline
%   ERR_PLT: If an error plot is requested (as the third output), then an error
%       plot is generated and its handle is returned.
%   NOTE -- If only 2 arguments are called for, then no plot is drawn.
%
%   TODO: For functions that produce more than 2 inputs, plots cannot be
%   visulaized. To accomodate this, we should take a set of 2D points and then
%   create higher-D points from them so that the values can be plotted.
%
% Example:
%   fx        = @(x) sin(2*pi*x);
%   bli_ip    = BLI(fx);
%   spline_ip = Spline(fx);
%   lg_ip     = Lagrange(fx);
%   n_tpts    = 1000;
%   bounds    = [-1; 1];
%   [err, sp] = compareIObjs(n_tpts, bounds, fx, bli_ip, 'BLI', spline_ip, ...
%       'Spline', lg_ip, 'Lagrange');
%
% Example: This can also be used to just analyze a single object
%   fx  = @(x) sin(2*pi*x);
%   n_tpts = 100;
%   bounds = [-1; 1];
%   compareIObjs(n_tpts, bounds, fx);
%
% Author: Archit Gupta, March 03, 2017
% Updated: October 24, 2017 (Updated documentation and added the functionality
%   for multiple interpolants to be compared at the same time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    multiple_interpolants = (nargin > 5); 
    no_interpolants       = (nargin == 3);

    if ( nargin < 5 )
        % Just comparing function against one interpolant and label not specified
        if (no_interpolants)
            ip_label = {};
        else
            ip_label = {'Interpolant'};
        end
    end

    % Check that an EVEN number of entries are present in VARARGIN
    if rem(length(varargin), 2)
        error('Expecting an even number of Interpolant/Label entries!')
    end

    n_interpolants = length(varargin)/2;
    ip_handles     = cell(n_interpolants, 1);
    ip_labels      = cell(n_interpolants, 1);

    f_handle       = getEvalHandle(f_obj);
    for ip = 1:n_interpolants
        ip_handles{ip} = getEvalHandle(varargin{ip*2-1});
        ip_labels{ip}  = varargin{ip*2};
    end


    n_in_dims = size(bounds, 2);
    n_op_dims = size(f_handle(zeros(n_in_dims, 1)), 1);

    % [TODO]: There is a bug in the code when we are asking for just 1 sample
    % points. Set n_test_pts = [1, 1] and the only point that you get is [NaN,
    % NaN].
    n_test_pts = 1 + randi(floor([n_tpts_seed ( n_tpts_seed + 11)]/n_in_dims), ...
        1, n_in_dims);

    g_bounds = zeros(2, n_in_dims);
    for d_i = 1:n_in_dims
        g_bounds(1, d_i) = bounds{d_i}(1);
        g_bounds(2, d_i) = bounds{d_i}(end);    % For piecewise interpolants
    end

    test_obj     = TestingPoints(n_in_dims, n_test_pts, scaleData(g_bounds, 0.99));
    n_test_cases = test_obj.getNPts();

    f_vals       = zeros([n_op_dims n_test_pts]);
    ip_vals      = cell(n_interpolants, 1);
    [ip_vals{:}] = deal(zeros([n_op_dims n_test_pts]));

    x_in         = cell(n_test_cases, 1);
    sub          = cell(n_test_cases, 1);
    % Getting all the points and sub-indices... This can be pretty
    % time-consuming

    for t_i = 1:n_test_cases
        [x_in{t_i}, sub{t_i}] = test_obj.getPtAt(t_i);
    end
    fprintf(2, 'Testing with %d sample points.\n', n_test_cases); 

    tic;
    for t_i = 1:n_test_cases
        f_vals(:, sub{t_i}{:})  = f_handle( x_in{t_i} );
        %Another way (somewhat more stressful) of looking at the timing information
        %f_anonymous             = @() f_handle(x_in);
        %f_eval_trials(sub{:}, 1) = timeit(f_anonymous);
    end
    f_eval_time = toc;
    fprintf(2, 'Function Evaluation time: %0.4f\n', f_eval_time);

    %{
    figure, plot(reshape(f_eval_times, [], 1), 'LineWidth', 2.0);
    title('Function Evaluation');
    xlabel('Iteration');
    ylabel('Time (in s)');
    grid on;
    set(gca, 'FontSize', 28);
    fprintf(2, 'Function Evaluation time: %0.8f\n', mean(reshape(f_eval_trials, [], 1)));
    %}
    
    ip_eval_times = zeros(n_interpolants, 1);
    for ip = 1:n_interpolants
        tic;
        for t_i = 1:n_test_cases
            ip_vals{ip}(:, sub{t_i}{:}) = ip_handles{ip}( x_in{t_i} );
            %c_anonymous               = @()  c_handle(x_in);
            %ip_eval_trials(sub{:}, 1) = timeit(c_anonymous);
        end
        ip_eval_times(ip) = toc;
        fprintf(2, '%s Evaluation time: %0.4f\n', ip_labels{ip}, ip_eval_times(ip));

        %{
        figure, plot(reshape(bli_eval_times, [], 1), 'LineWidth', 2.0);
        title('BLI Evaluation');
        xlabel('Iteration');
        ylabel('Time (in s)');
        grid on;
        set(gca, 'FontSize', 28);
        fprintf(2, '%s Evaluation time: %0.8f\n', mean(reshape(ip_eval_trials, [], 1)));
        %}
    end

    % Speedup calculation
    speedup = f_eval_time ./ ip_eval_times;

    % Error calculation
    est_err    = zeros(n_interpolants, 1);
    error_args = ErrorArgs();

    error_args.setDefaults();
    error_args.setMode('MAX_{REL}');
    error_args.setEps(1e-6);

    for ip = 1:n_interpolants
        est_err(ip) = find_error(ip_vals{ip}-f_vals, f_vals, error_args.eps(), error_args.errorMode());
        fprintf(2, '%s Error: %.10f\n', ip_labels{ip}, est_err(ip));
    end

    s_pts = cell(n_in_dims, 1);
    for d_i = 1:n_in_dims
        s_pts{d_i} = test_obj.getPts(d_i);
    end

    speedup = f_eval_time ./ ip_eval_times;

    if (nargout == 2)
        return;
    end

    fig_handle = plotComparison({}, s_pts, f_vals, ip_vals{:}, 'Baseline', ip_labels{:});

    if (nargout > 2)
        err_plt = plotDiffs({}, s_pts, f_vals, ip_vals{:}, ip_labels{:});
    end
end
