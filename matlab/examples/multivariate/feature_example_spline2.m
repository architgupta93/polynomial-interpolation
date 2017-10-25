function outs = feature_example_spline2(function_handle)
% function outs = feature_example_spline2(function_handle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use of Spline Interpolant for multi-variate function approximation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Example illustrating the use of Multi-variate Spline Interpolant.');

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
    args = defaultInterpolantArgs(2);

    disp('Instantiating 2D Spline.')
    bl_interpolant = Spline2D(function_handle, args{:});

    disp('Choosing random test points.')
    n_test_pts  = 1000;
    t_pts       = [sort(rand(1, n_test_pts)); sort(rand(1, n_test_pts))];
    t_vals      = function_handle(t_pts);

    interpolated_values = zeros(size(t_pts, 2));
    for t = 1 : n_test_pts
        interpolated_values(t) = bl_interpolant.computeWithDer(t_pts(:, t));
    end

    disp('Plotting the interpolated data against ground truth.')
    figure(); hold on;
    plot(t_pts(1,:), t_vals, 'LineWidth', 1.5);
    plot(t_pts(1,:), interpolated_values, 'LineWidth', 1.5)
    xlabel('test points');
    ylabel('values');
    legend('ground truth', 'interpolated');
    set(gca, 'FontSize', 28);
    grid on;
end

% Local functions here which are internal to the script
function f_handle = pickFunctionHandle();
    disp('Choosing function handle for testing.');
    f_handle = @(x) sin(2*pi*sqrt(x(1,:).^2 + x(2,:).^2));
end
