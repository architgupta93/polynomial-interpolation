function success = setupPPIPaths()
% function success = setupPPIPaths()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Helper script for automatically setting up the PATH variable to acces various
% scripts and classes defined in the package.
%
% If you would like to add a new interpolant class to the package, add all the
% contents in a flattened hierarchy in a folder, say  <foldername>, and add
% that folder to the main matlab directory containing this script. Its paths
% can be setup simply by adding <foldername> to the cell array
% 'path_locations_to_add'.
%
% OUTPUTS:
%   success: Boolean variable signifying that the paths were set up during this
%   call to the function
%
% TODO: Since we are introducing a persistent variable in this script, we
% should add a script to clear all the unwanted persistent/global variables
% created by the package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Archit Gupta
% Date: August 18, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name_of_this_file = 'setupPPIPaths.m';
    PPI_PATHS_SET_VALUE = 'PATHS_SET';
    persistent PPI_PATHS_SET; 

    if (strcmp(PPI_PATHS_SET, PPI_PATHS_SET_VALUE))
        disp('It looks like PPI paths have already been set. Aborting!')
        success = 0;
        return;
    else
        current_dir = strcat(pwd, '/');
        
        % Check that the function has been called from the directory in which
        % the .m file is saved. All paths will be used relative to the location
        % of this file.

        if (~exist(strcat(current_dir, name_of_this_file)))
            error('PATH setup script not called from its file location. Please CD to the correct location and try again.')
        else
            path_locations_to_add = {'common', ...
                'utils', ...
                'examples', ...
                'splines', ...
                'cosines', ...
                'cheb-series', ...
                'utils/test_functions', ...
                'bli'};
            for loc = 1 : length(path_locations_to_add)
                addpath(strcat(current_dir, path_locations_to_add{loc}));
            end
        end
    end
    PPI_PATHS_SET = PPI_PATHS_SET_VALUE;
    disp('PATHS have been set up succesfully for PPI.')
    success = 1;
end
