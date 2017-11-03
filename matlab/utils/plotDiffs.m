function success = plotDiffs(axis_labels, ref_pts_or_individual_refs, baseline, varargin)
% function success = plotDiffs(axis_labels, ref_pts_or_individual_refs, baseline, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the ERROR between a fixed 'baseline' and a variable number of reference
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
%   BASELINE: The baseline argument with respect to which, the error in other
%       measurements should be reported
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
%       plotDiffs({'time (s)', 'v(t)'}, {ts}, xs, ys, 'y');
%
%       Multi-Variate...
%       xs       = linspace(0, 1, 99);
%       ys       = linspace(0, 1, 101);
%       [YS, XS] = meshgrid(ys, xs);
%       f        = sin(2*pi*XS) + cos(2*pi*YS);
%       g        = cos(2*pi*XS) - sin(2*pi*YS);
%       plotDiffs({'x', 'y', 'v(x,y)'}, {xs, ys}, f, g, 'g');
%
%   See: plotComparison, compareIObjs, compareIDers
%
%   TODO: Add the functionality to treat the first data set as a baseline and
%   then interpolate the rest at the corresponding values (maybe, not sure how
%   this can be done)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Archit Gupta
% Date: November 02, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    success = 0;

    % Implementation for USAGE 01
    if rem(length(varargin), 2)
        error('Got %d arguments when expecting the same number of data entries as labels.', length(varargin));
    end

    n_comp_objs = length(varargin)/2;
    err_vals    = cell(n_comp_objs, 1);
    err_labels  = cell(n_comp_objs, 1);

    for ci = 1 : n_comp_objs
        try
            err_vals{ci}   = baseline - varargin{ci};
            err_labels{ci} = strcat('\Delta{', varargin{ci+n_comp_objs}, '}');
        catch ME
            % This is now how this should be done but I am in a hurry
            fprintf('DEBUG: %s\n', ME.message);
            error('Mismatch in array sizes to be diffed.. Come here, take a look!')
        end
    end

    % Again, USAGE 01 for plotComparison
    success = plotComparison(axis_labels, ref_pts_or_individual_refs, ...
        err_vals{:}, err_labels{:});
end
