disp('Using of Barycentric-Lagrange Interpolant for smoothing');
% Try this for a discontinuity
fun             = getTestFHandle(1, 1, 'd');

% This one has a derivative discontinuity
% fun             = getTestFHandle(1, 1, 'c');

disp('Generating arguments for the interpolant.');
args            = defaultInterpolantArgs();
args{4}         = 'chebyshev';

disp('Instantiating BLI.')
bl_interpolant  = Spline1D(fun, args{:});
% bl_interpolant.plotChebCoeffs();

hod_args        = args;
hod_args{3}(:)  = 6;
hod_interpolant = Spline1D(fun, hod_args{:});
% hod_interpolant.plotChebCoeffs();

n_test_pts      = 1000;
bounds          = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO');
