disp('Example illustrating the use of a Smoothing Interpolant.');
[fun, dfun]    = getTestFHandle(1, 1, 'd');

% step_width      = 1;
% fun             = @(x) abs(x) .* tanh(step_width*x);
% dfun            = @(x) (sign(x) .* tanh(step_width*x)) + ...
%                     (step_width * abs(x) .* (1 - tanh(step_width*x).^2));

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating SmoothStep.')
ss_interpolant = SmoothStep(fun, args{:});

% Get the derivatives
l_derv         = dfun(args{2}{1}(1));
r_derv         = dfun(args{2}{1}(end));
coeffs         = ss_interpolant.fit(l_derv, r_derv)

n_test_pts     = 1000;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, ss_interpolant);
