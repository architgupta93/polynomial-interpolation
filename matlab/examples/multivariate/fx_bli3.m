n_dims = 3;
        
disp('Example illustrating the use of Multi-variate Barycentric-Lagrange Interpolant.');
fun = getTestFHandle(n_dims, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args = defaultInterpolantArgs(n_dims);

disp('Instantiating 3D BLI.')
bl_interpolant = BLI3(fun, args{:});

n_test_pts = 100;
bounds     = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant);
