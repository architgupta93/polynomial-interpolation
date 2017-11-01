disp('SplineBLI interpolant for smoothing');
% Discontinuous function
[fun, dfun, bkp] = getTestFHandle(1,1,'d');

disp('Generating arguments for the interpolant.');
args            = defaultInterpolantArgs();
args{4}         = 'chebyshev';
disp('Instantiating SplineBLI.')
bl_interpolant  = SplineBLI(fun, args{:});
% bl_interpolant.plotChebCoeffs();

hod_args        = args;
hod_args{3}(:)  = 6;
hod_interpolant = SplineBLI(fun, hod_args{:});
% hod_interpolant.plotChebCoeffs();

n_test_pts      = 1000;
bounds          = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO');
