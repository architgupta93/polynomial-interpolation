disp('SplineBLI interpolant for smoothing');
% Discontinuous function
% [fun, dfun, bkp] = getTestFHandle(1,1,'d');

step_width      = 10000;
delta           = 0.0001;
fun             = @(x) (1 - x) .* (1 + tanh(step_width*(x-delta)))/2;
dfun            = @(x) -(1 + tanh(step_width*(x-delta)))/2 + ...
    (step_width/2)*(1-x).*(1 - tanh(step_width*(x-delta)).^2);

disp('Generating arguments for the interpolant.');
args            = defaultPiecewiseInterpolantArgs();
args{2}         = {[-1.0; -0.1; 0.1; 1.0]};
args{3}(:)      = 3;
args{4}         = 'chebyshev';
disp('Instantiating SplineBLI.')
% Pass in argument for smoothing
bl_interpolant  = PiecewiseBLI(fun, args{:}, false);

% Set the middle piece to be a low order spline
bl_interpolant.setInterpolant(2, @Spline1D, fun, {[-0.1; 0.1]}, [2]);

% bl_interpolant.plotChebCoeffs(1);
% bl_interpolant.plotChebCoeffs(2);

hod_args        = args;
hod_interpolant = PiecewiseBLI(fun, hod_args{:}, true);
% hod_interpolant.plotChebCoeffs(1);
% hod_interpolant.plotChebCoeffs(2);

n_test_pts      = 1000;
bounds          = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'PW-BLI', ...
    hod_interpolant, 'SPLINE-BLI');
