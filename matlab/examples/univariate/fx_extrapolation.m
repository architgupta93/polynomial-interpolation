disp('Example illustrating the Extrapolation with Splines.');
fun            = getTestFHandle(1, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating Passive Spline.')
ps_interpolant = PassiveSpline1D(fun, args{:});

disp('Instantiating Not-a-Knot Spline.')
de_interpolant = Spline1D(fun, args{:});

n_test_pts     = 1000;
bounds         = scaleData(args{2}, 1.5);
compareIObjs(n_test_pts, bounds, fun, ps_interpolant, 'PASSIVE', de_interpolant, 'SPLINE');