disp('Example illustrating the use of Barycentric-Lagrange Interpolant.');
fun            = getTestFHandle(1, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating BLI.')
bl_interpolant = BLI(fun, args{:});

n_test_pts     = 1000;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant);
