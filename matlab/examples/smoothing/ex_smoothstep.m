disp('Using of Barycentric-Lagrange Interpolant for smoothing');
step_width      = 100;
% [ofun, d_fun, bkpts] = getTestFHandle(1, 1, 'd'); 
% fun             = @(x) (1 - x) .* (1 + tanh(step_width*x))/2;
% dfun            = @(x) (step_width/2)*(1-x).*(1 - tanh(step_width*x).^2)  - (1 + tanh(step_width*x))/2;

fun             = @(x) (1 - x) .* (1 + tanh(step_width*x))/2;
dfun            = @(x) (step_width/2)*(1-x).*(1 - tanh(step_width*x).^2)  - (1 + tanh(step_width*x))/2;

disp('Generating arguments for the interpolant.');
args            = defaultInterpolantArgs();
args{4}         = 'chebyshev';

disp('Instantiating BLI.')
bl_interpolant  = BLI(fun, args{:});
% bl_interpolant.plotChebCoeffs();

hod_args        = args;
hod_args{3}(:)  = 6;
hod_interpolant = BLI(fun, hod_args{:});
% hod_interpolant.plotChebCoeffs();

pw_args         = defaultPiecewiseInterpolantArgs(1);
pw_args{3}(:)   = 5;
pw_args{4}      = 'chebyshev';
pw_interpolant  = PiecewiseBLI(fun, pw_args{:});

n_test_pts      = 1000;
bounds          = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO', pw_interpolant, 'PW-BLI');
compareIDers(n_test_pts, bounds, dfun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO', pw_interpolant, 'PW-BLI');

