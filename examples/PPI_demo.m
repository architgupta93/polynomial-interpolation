function outs = illustrative_example(function_handle)
% function outs = illustrative_example(function_handle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ILLUSTRATIVE EXAMPLE on the use of (piecewise) polynomial interpolants
% TODO: Make the function interactive for easy first time use
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Illustrative example on the use of polynomial interpolation.');

    % First, we check if the user supplied a function handle
    if (nargin > 0)
        if (isa(function_handle, 'function_handle'))
            disp('User supplied function handle %s', function_handle);
        else
            disp('Invalid input supplied, expecting function handle.');
            function_handle = pickFunctionHandle();
        end
    else
        function_handle = pickFunctionHandle();
    end

    disp('Generating arguments for the interpolant.');
    args = defaultInterpolantArgs();

    disp('Instantiating a SPLINE Interpolant.')
    spline_interpolant = Spline1D(function_handle, args{:});

    disp('Choosing random test points.')
    n_test_pts  = 100;
    t_pts       = sort(rand(n_test_pts, 1));
    t_vals      = function_handle(t_pts);

    interpolated_values = zeros(size(t_pts));
    for t = 1 : length(t_pts)
        interpolated_values(t) = spline_interpolant.computeWithDer(t_pts(t));
    end

    disp('Plotting the interpolated data against ground truth.')
    figure(); hold on;
    plot(t_pts, t_vals, 'LineWidth', 1.5);
    plot(t_pts, interpolated_values, 'LineWidth', 1.5)
    xlabel('test points');
    ylabel('values');
    legend('ground truth', 'interpolated');
    set(gca, 'FontSize', 28);
    grid on;
end

% Local functions here which are internal to the script
function f_handle = pickFunctionHandle();
    disp('Choosing function handle for testing.');
    f_handle = @(x) sin(2*pi*x);
end
