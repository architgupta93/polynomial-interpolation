disp('Example illustrating the use of Barycentric-Lagrange Interpolant.');
fun = getTestFHandle(1, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args = defaultInterpolantArgs();

disp('Instantiating Lagrange Interpolant.')
lag_interpolant = Lagrange(fun, args{:});

n_test_pts = 100;
bounds     = args{2};
compareIObjs(n_test_pts, bounds, fun, lag_interpolant);
