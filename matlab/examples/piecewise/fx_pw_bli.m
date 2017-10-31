n_dims         = 1;
n_pieces       = 4;

disp('Example illustrating the use of Barycentric-Lagrange Interpolant.');
fun            = getTestFHandle(1, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args           = defaultPiecewiseInterpolantArgs(1);

disp('Instantiating BLI.')
bl_interpolant = PiecewiseBLI(fun, args{:});

n_test_pts     = 100;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant);
