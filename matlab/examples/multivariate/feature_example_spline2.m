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
    sp_interpolant = Spline2D(function_handle, args{:});

    n_test_pts = 1000;
    bounds     = {[-1; 1], [-1; 1]};
    compareIObjs(n_test_pts, bounds, function_handle, sp_interpolant);
end

% Local functions here which are internal to the script
function f_handle = pickFunctionHandle();
    disp('Choosing function handle for testing.');
    f_handle = @(x) sin(2*pi*sqrt(x(1,:).^2 + x(2,:).^2));
end
