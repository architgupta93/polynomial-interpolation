disp('Example illustrating the use of Barycentric-Lagrange Interpolant.');
% fun            = getTestFHandle(1, 1, 'smooth');
step_width      = 10000;
fun             = @(x) tanh(step_width*x);
dfun            = @(x) step_width*(1 - tanh(step_width*x).^2);

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating BLI.')
bl_interpolant = SmoothStep(fun, args{:});

n_test_pts     = 1000;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant);
