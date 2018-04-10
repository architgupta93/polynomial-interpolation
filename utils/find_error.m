function [calculated_error, errorPts] = find_error(err, sol, threshold, ...
    MODE, doPlots)
% function [calculated_error, errorPts] = find_error(err, sol, threshold, ...
%     MODE, doPlots)
% Author: Archit Gupta, April 22, 2016
% Modified (November 11, 2016): Added plotting functionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script evaluates the "ERROR" metric for a given data. Several
% metrics/methods can be specified to calculate the error.
% INPUTS:
%   err -   A matrix/vector of error values
%   sol -   The actual function values (the deviation from which has been
%           reported as the error    
%   threshold - any value which is less than threshold*max(sol) will be ignored
%   for computing the error
%           that should be ignored while calculating the error
%   MODE -  What error metric is required (MAX/MEAN/MEDIAN) - (ABS/REL). Some
%           combinations are allowed at present
%   doPlots - (Optional) Visualizing the points where the error is significant
%
% OUTPUTS:
%   calculated_error - Calculated value of the ERROR metric
%   errorPts        - Currently just gives the 1 point at which the error metric
%                     has the maximum value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin < 5)
        doPlots = 0;
    end

    calculated_error = Inf;
    % Get the dimensions of the solution matrix
    sol_dims = size(sol);
    errorPts = Inf*ones(size(sol_dims));

    err_val = reshape(abs(err), [], 1);
    sol_val = reshape(abs(sol), [], 1);
    max_sol_val = max(sol_val);
    eps = threshold*max_sol_val;    % Anything that is 4 orders of magnitude less
                                    % than the max in ignored 
    if (length(err_val) ~= length(sol_val))
        fprintf(2, 'Length of Error and Solution matrices is not consistent\n');
        return;
    end

    sol_val = sol_val.*(sol_val > eps);
    err_val = err_val.*(sol_val > eps); % if sol_val ~ 0, ignore the error
    error_pts_base_val = Inf;  % Since we have reshaped the original tensor into a
                            % 1D array, some extra effort will be required to
                            % get the original corrdinates (some remainders and
                            % substractions will be required). error_pts_base_val
                            % is the index in the RESHAPED 1D ARRAY 
    pointwise_error = err_val./sol_val;
    % Default error provided by the script is MAXIMUM RELATIVE ERROR
    [calculated_error, error_pts_base_val] = nanmax(pointwise_error);
    if (strcmp(MODE, 'MAX_{REL}'))
        % NOTHING TO DO: This is the default implementation
    elseif (strcmp(MODE, 'MAX_{ABS}'))
        [calculated_error, error_pts_base_val] = nanmax(err_val);
    elseif (strcmp(MODE, 'MEAN_{REL}'))
        calculated_error = nanmean(pointwise_error);
    elseif (strcmp(MODE, 'MEAN_{ABS}'))
        calculated_error = nanmean(err_val);
        [~, error_pts_base_val] = nanmax(err_val);
    elseif (strcmp(MODE, 'MEDIAN_{REL}'))
        calculated_error = nanmedian(pointwise_error);
    elseif (strcmp(MODE, 'MEDIAN_{ABS}'))
        calculated_error = median(err_val);
        [~, error_pts_base_val] = nanmax(err_val);
    else
        fprintf(2, 'ERROR: Error MODE Unknown!\n');
        return;
    end 

    % There might be some skew in the data (one/few very large values, which
    % make us omit pretty much all the values in the solution matrix while
    % computing the error. Let's print out the number of entries that were
    % considered for giving out the error metric

    if (length(strfind(MODE,'REL')) ~= 0)
        % i.e. If you were looking for a relative error metric
        nans_in_relative_error = isnan(pointwise_error);
        nNaNs = sum(nans_in_relative_error);
        fprintf(2, 'Calculated Relative Error ignoring %d of %d points\n', ...
            nNaNs, length(sol_val));
    end 

    if (nargout > 1)    % No Extra work if not required
        for i = 1:length(sol_dims)  % Get the index along each dimension
            errorPts(i) = floor(error_pts_base_val/sol_dims(i));
            error_pts_base_val = error_pts_base_val - sol_dims(i)*errorPts(i);
        end
    end

     % Show a plot of the points that have a relatively high error. For now, we
     % will just hardcode this error to be 5e-2, i.e., 5% error -- Makes sense
     % for relative error calculation only 
     HIGH_STAKE_ERROR_THRESH = 5e-2;

    if (doPlots)   
        high_stake_error_points = nan(sol_dims);    % NaNs don't show up plots
        high_stake_error_points(:) = (pointwise_error>HIGH_STAKE_ERROR_THRESH);     
        plot_line_width = 2.5;
        if (length(sol_dims) == 1)
            figure(), plot(sol, 'LineWidth', plot_line_width);
            hold on;
            % Scatter with appropriate scaling so that both the plots are
            % visible (Just multplying with the max value for now
            stem(high_stake_error_points*max_sol_val, 'linestyle', 'none', ...
                'LineWidth', plot_line_width);
            set(gca, 'FontSize', 28);
            legend('Original Solution', 'High error points');
        elseif (length(sol_dims) == 2)
            figure(), surf(reshape(sol, sol_dims));
            freezeColors;
            hold on;
            stem3(high_stake_error_points*max_sol_val, 'LineWidth', ...
                plot_line_width);
            set(gca, 'FontSize', 28);
            legend('Original Solution', 'High error points');
        else
            fprintf(['ERROR: Too many dimensions (%d) in the input data for', ...
                'visualization!\n'], length(sol_dims));
        end
    end
end
