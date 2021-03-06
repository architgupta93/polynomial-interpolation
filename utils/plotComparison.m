function success = plotComparison(axis_labels, ref_pts_or_individual_refs, varargin)
% function success = plotComparison(axis_labels, ref_pts_or_individual_refs, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot multiple overlaid line or surf plots for comparison.
% INPUT(s):
%   AXIS_LABELS: List of axis labels in the form of a cell array.
%       - (x,y) for 1D plots, and
%       - (x,y,z) for 2D plots.
%       This could optionally be passed in as {}, in which case some defaults
%       would be chosen
%
%   REF_PTS_OR_INDVIDUAL_REFS: Two inputs are possible:
%       1. Reference points, i.e., points on the axis to which the data
%       corresponds. This could also be left empty and something appropriate
%       will be chosen, however, in this case, it will be assumed that all
%       subsequent entries are DATA, LABEL entries, and have been sampled at
%       the same points.
%       
%       2. [], i.e., individual DATA, LABEL entries will also be accompanied
%       with the points at which they were sampled. In this all values will be
%       interpolated at the baseline ref points and plotted
%
%   VARARGIN: variable argument list (must have either of the following
%       arrangement, depending on REF_PTS_OR_INDVIDUAL_REFS specification
%       above): 
%       {<data1>, ..., <datan>, <label1>, ..., <labeln>}, OR
%       {<ref_pts1>, ..., <ref_ptsn>, <data1>, ..., <datan>, <label1>, ..., <labeln>}
%
%   Example USAGE 01:
%       ts = linspace(0, 1, 100);
%       xs = sin(2*pi*ts);
%       ys = cos(2*pi*ts);
%       plotComparison({'time (s)', 'v(t)'}, {ts}, xs, ys, 'x', 'y');
%
%       Multi-Variate...
%       xs       = linspace(0, 1, 99);
%       ys       = linspace(0, 1, 101);
%       [YS, XS] = meshgrid(ys, xs);
%       f        = sin(2*pi*XS) + cos(2*pi*YS);
%       g        = cos(2*pi*XS) - sin(2*pi*YS);
%       plotComparison({'x', 'y', 'v(x,y)'}, {xs, ys}, f, g, 'f', 'g');
%
%   Example USAGE 02:
%       xt = linspace(0, 1, 100);
%       xs = sin(2*pi*x);
%       yt = linspace(0, 2, 101);
%       ys = cos(2*pi*x);
%       plotComparison({'time (s)', 'v(t)'}, [], xt, yt, xs, ys, 'x', 'y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    success = 0; % Look at the pessimism

    if (isempty(ref_pts_or_individual_refs))
        % Usage 2
        error('Usage 2 in the documentation has not been implemented yet!')

        % TODO: When implemented, keep the return statement here.
        return;
    end

    ref_pts = ref_pts_or_individual_refs;

    if rem(length(varargin), 2)
        error('Got %d arguments when expecting the same number of data entries as labels.', length(varargin));
    end

    n_overlays  = length(varargin)/2;
    % Try to find the nature of plotting to be used...
    % PLOT for 1D,
    % SURF/MESH for 2D
    data_labels = {varargin{n_overlays+1:end}};
    data_vals   = cell(n_overlays, 1);

    for ol = 1:n_overlays
        data_vals{ol} = varargin{ol};
    end

    % Get the number of data dimensions to be plotted as well as the number of plots required.
    n_in_dims    = length(ref_pts);
    is_1D_data   = (n_in_dims == 1);

    % These values have to be sorted otherwise, we get a dimension mismatch in
    % plotting. Also, we do not want repeatitions, hence randperm
    var_in_dims  = sort(randperm(n_in_dims, min(2, n_in_dims)), 'ascend');

    in_cols      = cell(1, n_in_dims);
    [in_cols{:}] = deal(':');

    out_dims = ndims(data_vals{1}) - n_in_dims;
    out_cols = cell(1, out_dims);

    % Pick a random dimension to be plotted.
    % We already assume that the first out_dims dimensions correspond to the outputs
    data_size = size(data_vals{1});
    for o_dim = 1:out_dims
        if data_size(o_dim) > 1
        out_cols{o_dim} = randi([1 data_size(o_dim)]);
        else
            out_cols{o_dim} = 1;
        end
    end

    % Print the output dimension that will be plotted.
    fprintf(2, 'Plotting Output dimension: ');
    disp([out_cols{:}]);

    if (~is_1D_data)

        if (g_dims(data_vals{1} > 2))
            % If this is the case, then we can only plot a subset of the indices
            % that are plotted. The current scheme is simple!
            % Just pick all but 2 dimensions at random and fix indices for all
            % the remaning dimensions.
            
            for in_dim = 1 : n_in_dims
                if ( (in_dim == var_in_dims(1)) || (in_dim == var_in_dims(2)) ) 
                    % Let all elements be allowed to vary for these dimensions
                    continue;
                end

                in_cols{in_dim} = randi([1 length(ref_pts{in_dim})], 1);
            end
        end

        % TODO: What we did for input colons can also be done for output colons
        % but for the time being we will just stick to something simple
        if (out_dims > 0)
            % DO SOMETHING
        end
    end

    % Even if the data is NOT 1D, we can simply plot the indices of the sample
    % points to look at the error (instead of actually plotting the
    % sample-point, sample-value pairs). However, that is NOT what is done right
    % now: TODO
    is_1D_plot = is_1D_data;

    figure(); hold on;
    if (is_1D_plot)
        % Just change these and the rest should work out
        % For example:
        %   plot_fun  = @scatter;
        %   plot_args = {};
        % Should produce a scatter plot instead
        plot_fun  = @plot;
        plot_args = {'LineWidth', 2.0};

        if (isempty(axis_labels))
            axis_labels = {'pt', 'v(pt)'};
        end

        % Loop thorugh the data and plot
        % TODO: Maybe put a try catch block around this.
        for ol = 1:n_overlays
            plot_fun(ref_pts{1}, squeeze(data_vals{ol}(out_cols{:}, in_cols{:})), plot_args{:});
        end
        xlabel(axis_labels{1});
        ylabel(axis_labels{2});
    else
        plot_fun  = @surf;
        plot_args = {};

        if (isempty(axis_labels))
            axis_labels = {'pt1', 'pt2', 'v(pt1, pt2)'};
        end

        for ol = 1:n_overlays
            % TODO: The requirement of a transposing the data for a surf plot
            surf(ref_pts{var_in_dims(1)}, ref_pts{var_in_dims(2)}, squeeze(data_vals{ol}(out_cols{:}, in_cols{:}))');
        end
        xlabel(axis_labels{var_in_dims(1)});
        ylabel(axis_labels{var_in_dims(2)});
    end
    legend(data_labels{:});
    set(gca, 'FontSize', 28);
    grid on;
    hold off;

    % TODO: What about higher dimensions? For ever 1/2, can we let the user
    % supply arguments and handles for plotting
    success = 1;
end
