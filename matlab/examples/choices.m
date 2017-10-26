% Here is a simple experiment to decide how multi-dimensional matrices should be used.
% Let us start with a fairly large matrix, and time some anonymous functions
t_mat = rand(5, 7, 9, 11, 10, 17, 2, 3, 5);
f1 = @() t_mat(1, :, :, :, :, :, :, :, :) - t_mat(2, :, :, :, :, :, :, :, :);
f2 = @() t_mat(:, :, :, :, :, :, :, :, 1) - t_mat(:, :, :, :, :, :, :, :, 2);

% Now try:
timeit(f1)
timeit(f2)

% This will give a good idea of which of the two methods is more appropriate
% f2 should take significantly lesser time that f1.
% This is because of the way multi-dimensional matrices are stored in MATLAB
