n_dims         = 3;
n_pieces       = 4;

disp('Example illustrating the use of Barycentric-Lagrange Interpolant.');
fun            = getTestFHandle(n_dims, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args           = defaultPiecewiseInterpolantArgs(n_dims);

disp('Instantiating BLI.')
bl_interpolant = PiecewiseBLI3(fun, args{:});

n_test_pts     = 100;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'Piecewise BLI');