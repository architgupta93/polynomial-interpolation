%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ILLUSTRATIVE EXAMPLE on the use of (piecewise) polynomial interpolants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outs = illustrative_example(function_handle)

    disp('Illustrative example on the use of polynomial interpolation.');

    % First, we check if the user supplied a function handle
    if (nargin > 0)
        if (isa(function_handle, 'function_handle'))
            disp('User supplied function handle %s', function_handle);
        else
            disp('Invalid input supplied, expecting function handle.');
            function_handle = pickFunctionHandle();
        end
    else
        function_handle = pickFunctionHandle();
    end

    disp('')
end

% Local functions here which are internal to the script
function f_handle = pickFunctionHandle();
    disp('Choosing function handle for testing.');
    f_handle = @(x) sin(x);
end
