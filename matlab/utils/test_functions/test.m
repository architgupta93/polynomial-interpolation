clc; close all;
n_elems = 2;
n_test_pts = [50; 51];  % Points along each dimension
c_vec = 4*rand(n_elems, 1);
w_vec = 1*rand(n_elems, 1);

% f_handle = f_continuous(n_elems, 1, [w_vec; c_vec]);
% f_handle = f_discontinuous(n_elems, 1, [w_vec(1:2); c_vec]);
% f_handle = f_gaussian(n_elems, 1, [w_vec; c_vec]);
% f_handle = f_oscillatory(n_elems, 1, [w_vec(1); c_vec]);
f_handle = f_product_peak(n_elems, 1, [w_vec; c_vec]);
% f_handle = f_corner_peak(n_elems, 1, c_vec);

f_vals = zeros(n_test_pts');        % Need row vector here
t_args = TestingArgs(n_elems, n_test_pts);

for pt = 1:prod(n_test_pts)
    t_pt = t_args.getTestPoint(pt);
    f_vals(pt) = f_handle.evaluate(t_pt);
end

t_pts = cell(n_elems, 1);
for e = 1:n_elems 
    t_pts{e} = t_args.getSamplePoints(e);
end

% Plotting
figure(); hold on;
if (n_elems == 1)
    plot(t_pts{1}, f_vals, 'LineWidth', 2);
elseif (n_elems == 2)
    surf(t_pts{1}, t_pts{2}, f_vals');
    colormap autumn;
else
    fprintf('Cannot plot 3-or-higher dimensional data!\n');
end
set(gca, 'FontSize', 28);
grid on;
