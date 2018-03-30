disp('Example illustrating the use of Multi-variate Barycentric-Lagrange Interpolant.');
fun = getTestFHandle(2, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args = defaultInterpolantArgs(2);

disp('Instantiating 2D BLI.')
bl_interpolant = BLI2(fun, args{:});

n_test_pts = 100;
bounds     = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'BLI 2D');