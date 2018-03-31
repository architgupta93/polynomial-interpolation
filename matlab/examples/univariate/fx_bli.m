disp('Example illustrating the use of Splines and Barycentric-Lagrange Interpolant.');
fun            = getTestFHandle(1, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating BLI.')
bl_interpolant = BLI(fun, args{:});

disp('Instantiating Spline.')
sp_interpolant = Spline1D(fun, args{:});

n_test_pts     = 1000;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'BLI', sp_interpolant, 'SPLINE');