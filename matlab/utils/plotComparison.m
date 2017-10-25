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
%       plotComparison({'time (s)', 'v(t)'}, ts, xs, ys, 'x', 'y')
%
%   Example USAGE 02:
%       xt = linspace(0, 1, 100);
%       xs = sin(2*pi*x);
%       yt = linspace(0, 2, 101);
%       ys = cos(2*pi*x);
%       plotComparison({'time (s)', 'v(t)'}, [], xt, yt, xs, ys, 'x', 'y')
%
%   TODO: Add the functionality to treat the first data set as a baseline and
%   then interpolate the rest at the corresponding values (maybe, not sure how
%   this can be done)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    n_overlays = length(varargin)/2;
    % Try to find the nature of plotting to be used...
    % PLOT for 1D,
    % SURF/MESH for 2D

    figure(); hold on;
    if (g_dims(varargin{1}) == 1)
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
            plot_fun(ref_pts, varargin{ol}, plot_args{:});
        end
    elseif (g_dims(varargin{2} == 2))

        if (isempty(axis_labels))
            axis_labels = {'pt1', 'pt2', 'v(pt1, pt2)'};
        end

        error('Two dimensions not implemented yet!');
    else
        error('Too many data dimensions to plot!');
    end
    legend(varargin{n_overlays+1:end});
    set(gca, 'FontSize', 28);
    xlabel(axis_labels{1});
    ylabel(axis_labels{2});
    grid on;
    hold off;

    % TODO: What about higher dimensions? For ever 1/2, can we let the user
    % supply arguments and handles for plotting
end
