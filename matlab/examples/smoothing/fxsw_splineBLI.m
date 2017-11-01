disp('SplineBLI interpolant for smoothing');
% Discontinuous function
% [fun, dfun, bkp] = getTestFHandle(1,1,'d');

step_width      = 10000;
fun             = @(x) (1 - x) .* (1 + tanh(step_width*x))/2;
dfun            = @(x) (step_width/2)*(1-x).*(1 - tanh(step_width*x).^2)  - (1 + tanh(step_width*x))/2;

disp('Generating arguments for the interpolant.');
args            = defaultPiecewiseInterpolantArgs();
args{2}{1}
args{4}         = 'chebyshev';
smoothing       = false;
disp('Instantiating SplineBLI.')
% Pass in argument for smoothing
bl_interpolant  = PiecewiseBLI(fun, args{:}, smoothing);
% bl_interpolant.plotChebCoeffs();

hod_args        = args;
hod_args{3}(:)  = 6;
hod_interpolant = PiecewiseBLI(fun, hod_args{:}, smoothing);
% hod_interpolant.plotChebCoeffs();

n_test_pts      = 1000;
bounds          = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'LO', ...
    hod_interpolant, 'HO');
