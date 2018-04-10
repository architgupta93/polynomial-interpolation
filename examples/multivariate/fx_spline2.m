disp('Example illustrating the use of Multi-variate Spline Interpolant.');
fun = getTestFHandle(2, 2, 'smooth');

disp('Generating arguments for the interpolant.');
args = defaultInterpolantArgs(2);

disp('Instantiating 2D Spline.')
sp_interpolant = Spline2D(fun, args{:});

n_test_pts = 100;
bounds     = args{2};
compareIObjs(n_test_pts, bounds, fun, sp_interpolant, 'Spline 2D');