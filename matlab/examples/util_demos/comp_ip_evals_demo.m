% Sample script for demonstrating the Compare-Interpolant-Evaluation utility

% Define a function to work with
fx     = @(x) sin(2*pi*x);
df_dx  = @(x) 2*pi*cos(2*pi*x);

% Comparing function evaluation against a single interpolant
args   = defaultInterpolantArgs(1);
bli_ip = BLI(fx, args{:});
bli_lb = 'BLI';

n_tpts = 100;
bounds = {[-1; 1]};
compareIObjs(n_tpts, bounds, df_dx, bli_ip, bli_lb);
