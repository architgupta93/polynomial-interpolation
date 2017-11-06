disp('Example illustrating the use of a 2-dimensional Smoothing Interpolant.');
[fun, dfun] = getTestFHandle(2, 1, 'd');

step_width  = 1;
fun         = @(x) abs(x(1,:).*x(2,:)) .* tanh(step_width*(x(1,:).*x(2,:)));
dfun        = @(x) [(sign(x(1,:).*x(2,:)) .* x(2,:) .* tanh(step_width*(x(1,:).*x(2,:)))) + ...
                (step_width * x(2,:) .* abs(x(1,:).*x(2,:)) .* (1 - tanh(step_width*(x(1,:).*x(2,:))).^2)); ...
            (sign(x(1,:).*x(2,:)) .* x(1,:) .* tanh(step_width*(x(1,:).*x(2,:)))) + (step_width * ...
                x(1,:) .* abs(x(1,:).*x(2,:)) .* (1 - tanh(step_width*(x(1,:).*x(2,:))).^2))];

disp('Generating arguments for the interpolant.');
args        = defaultInterpolantArgs(2);

disp('Instantiating 2D SmoothStep.');
ss_itp      = SmoothStep2(fun, args{:});

% TODO: Get the derivatives and fit them

n_tpts      = 100;
bounds      = args{2};
compareIObjs(n_tpts, bounds, fun, ss_itp);
