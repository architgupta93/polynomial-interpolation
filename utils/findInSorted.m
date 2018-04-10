function right_bound = findInSorted(x, range)
% This script has been borrowed from
% http://stackoverflow.com/questions/20166847/faster-version-of-find-for-sorted-vectors-matlab
% and modified for our use. It takes in a RANGE which has been sorted in
% DESCENDING ORDER and a 'SINGLE VALUE' x, that we want to look for in the
% RANGE
    
    if ( ~isscalar(x) )
        error('ERROR: Expecting a scalar input. Aborting');
    end
    
    % Expect range to be a 1D array, which has been sorted in a descending
    % order. Its just specific to this package. Sorry if it has too many
    % restrictions to be used elsewhere
    
    lower_bound = range(end);
    upper_bound = range(1);
    
    % The bound indices have been named with reference to the real number
    % line
    left_bound = numel(range);
    right_bound = 0;
    loc = 1;  
    
    if (lower_bound >= x)
       loc = left_bound;
       right_bound = left_bound;
    elseif (upper_bound <= x)
        loc = right_bound;
        left_bound = right_bound;
    end
    
    while ( (left_bound-1 > loc) || (right_bound+1 < loc) )
        loc = (floor((right_bound + left_bound)/2));
        if (x < range(loc))
            right_bound = loc;
        else
            left_bound = loc;
        end
    end
end
