disp('Using of Barycentric-Lagrange Interpolant for smoothing');
step_width      = 10000;
% [fun, dfun, bkpts] = getTestFHandle(1, 1, 'd'); 
% fun             = @(x) tanh(step_width*x);
% dfun            = @(x) step_width*(1 - tanh(step_width*x).^2);

fun             = @(x) (1 - x) .* (1 + tanh(step_width*x))/2;
dfun            = @(x) (step_width/2)*(1-x).*(1 - tanh(step_width*x).^2)  - (1 + tanh(step_width*x))/2;

disp('Generating arguments for the interpolant.');
args            = defaultInterpolantArgs();
args{4}         = 'chebyshev';

disp('Instantiating BLI.')
bl_interpolant  = BLI(fun, args{:});
bl_interpolant.plotChebCoeffs();

hod_args        = args;
hod_args{3}(:)  = 6;
hod_interpolant = BLI(fun, hod_args{:});
hod_interpolant.plotChebCoeffs();

pw_args         = defaultPiecewiseInterpolantArgs(1);
pw_args{2}      = {[-1.0; -0.1; 0.0; 0.1; 1.0]};
pw_args{3}      = [2];
pw_args{4}      = 'chebyshev';
pw_bli          = PiecewiseBLI(fun, pw_args{:});
pw_bli.plotChebCoeffs(2);

pw_spline       = PiecewiseSpline(fun, pw_args{:});

n_test_pts      = 1000;
bounds          = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO', pw_bli, 'PW-BLI', ...
    pw_spline, 'PW-SPLINE');
compareIDers(n_test_pts, bounds, dfun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO', pw_bli, 'PW-BLI', ...
    pw_spline, 'PW-SPLINE');

