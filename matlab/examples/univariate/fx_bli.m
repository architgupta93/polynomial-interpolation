disp('Example illustrating the use of Splines and Barycentric-Lagrange Interpolant.');
% We create a vector function. It takes 1 input and produces 3 outputs. The
% interpolant can handle this easily without any input from the user.
fun            = getTestFHandle(1, 5, 'smooth');
%                               ^  ^
%                           Input  Output
%                      Dimensions  Dimensions

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating BLI.')
bl_interpolant = BLI(fun, args{:});

disp('Instantiating Spline.')
sp_interpolant = Spline1D(fun, args{:});

n_test_pts     = 1000;
bounds         = args{2};
compareIObjs(n_test_pts, bounds, fun, bl_interpolant, 'BLI', sp_interpolant, 'SPLINE');