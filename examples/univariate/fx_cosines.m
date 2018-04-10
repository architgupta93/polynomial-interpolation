disp('Example illustrating the use of DCT as an Interpolant.');
fun = getTestFHandle(1, 1, 'smooth');

disp('Generating arguments for the interpolant.');
args = defaultInterpolantArgs();

% Set the order of the current interpolation is not good enough
args{3} = 6;

disp('Instantiating Cosine interpolant.')
dct_interpolant = DCTI(fun, args{:});

n_test_pts = 100;
bounds     = args{2};
compareIObjs(n_test_pts, bounds, fun, dct_interpolant, 'DCT');