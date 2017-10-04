classdef Interpolant < SaveLoad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classdef Interpolant < SaveLoad
% This class defines an abstract "Interpolant" object and defines all the
% access functions that this class of objects should support
%
% Abstract class, CANNOT BE INSTANTIATED

% All the interpolants that have been implemented in this codebase inherit this
% class and have to support the API functions defined here to be instantiable.
%
% See also: PiecewiseInterpolant, BLI, Spline, Lagrange, DCTI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        % INPUT: As of now, we only allow functions of vectors to be
        % interpolated, a scalar number, therefore, defines its size
        in_dims = 0;

        % OUTPUT: This can have an arbitrary dimension, as a result, the
        % dimensionality needs to be stored as a row vector
        op_dims = 0;

        % REFINENESS: For each interpolant, there is a notion of 'refinement'.
        % ORDER defines the degree of refineness. Higher the order, larger the
        % size of the interpolant, and hopefully, better the interpolation. The
        % exact nature of this depends on the interpolant
        order = [];

        % WEIGHTS: All the interpolants are 'LINEAR', in the sense that,
        % scaling the underlying function's sample values scales the
        % interpolant and adding two functions leads to addition of the
        % interpolant. So, it can be defined in terms of some weights.
        wts = [];

        % SAMPLE POINTS: Points at which the sample values were provided.
        i_pts = InterpolationPoints();
        
        % DOMAIN: Domain bounds inside which the interpolant CAN interpolate the values
        bounds = cell(0);

        % COLONS: This is a MATLAB/OCTAVE centric field required for indexing
        % all the OUTPUTs corresponding to a specific (partial) input value.
        % For example, if the output has dimesions 2x3x4x5, then for any input
        % value, these values are indexed using {:, :, :, :}, which is stored
        % in colons
        colons = cell(0);

        % EXTRAPOLATION SLOPE: Outside the sampled domain, the extrapolation
        % can be specified either automatically or by the user by specifying
        % the extrapolation slope
        extrap_slope = [];

        % SAMPLE VALUES: function values that were sampled for constructing the
        % interpolant. Sometimes, these come in handy for interpolation.
        f_vals = [];
    end

    methods (Static = true, Access = protected)
        function success = isColumnVector(x)
        % Internal function
        % Functions here can only take in vector inputs. Specifically, the
        % input should be a column vector and this function is used to verify
        % that it is indeed the case.
            n_dims_in_x = length(size(x));
            n_cols_in_x  = size(x, 2);
            if ( (n_cols_in_x == 1) && (n_dims_in_x == 2) )
                success = 1;
            else
                success = 0;
            end
        end
    end

    methods (Access = public) % Abstract methods from MATLAB
        function [val, der] = computeWithDer(Obj, x_in, coeffs);
        % Primary function called for interpolation
        % Provides the interpolated value and the derivative at the point
            var = [];
            der = [];
        end
        function der = secondDerivative(Obj, x_in);
        % Provides the second derivative at the query point
            val = [];
            der = [];
        end
    end

    methods (Access = public)
        function Obj = Interpolant(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % Class Constructor
        %
        % Only called through child classes. Do not call directly
        % INPUTS:
        %   f_vals: Function values at a set of points, or, function_handle
        %       which can be evaluated.
        %   in_dims: Dimensionality of the input vector (SCALAR)
        %   bounds: Domain over which the interpolant is supposed to operate.
        %       Currently, only rectangular boundaries are supported. This is
        %       expected to be a CELL ARRAY of size in_dims x 1. Each entry is
        %       supposed to have 2 values (start, end).
        %   order: READ order description for class variables. This is
        %       interpolant specific.
        %   i_type_or_x_vals: Specify the TYPE of interpolation points or the
        %       actual value of the sample points, at which it should (or has
        %       been) be evaluated.
        %       If f_vals is a set of values, then we additionally need to
        %       match that the number of values and number of points is correct.
        %       Some of these features are [TODO].

            if (nargin == 0)
                % DEBUG MESSAGE [TODO]: Remove this later on
                % fprintf(2, 'Instantiating an EMPTY interpolant of class: %s! \n', class(Obj)); 
                return;
            end

            % Checking if interpolation type (i_type) has been provided
            if ( isstr(i_type_or_x_vals) )
                if ( strcmp(i_type_or_x_vals, 'uniform') || ...
                    strcmp(i_type_or_x_vals, 'u') )
                    Obj.i_pts = UniformPoints(in_dims, order, bounds);
                elseif ( strcmp(i_type_or_x_vals, 'chebyshev') || ...
                    strcmp(i_type_or_x_vals, 'c') )
                    Obj.i_pts = ChebyshevPoints1Analyst(in_dims, order, bounds);
                elseif ( iscell(i_type_or_x_vals) )
                    Obj.i_pts = InterpolationPoints(1, order, bounds, ...
                        i_type_or_x_vals)
                    Obj.wts = wts;

                else
                    error('ERROR: Supplied "type": %s Invalid!\n', ...
                        i_type_or_x_vals);
                end
            else
                error(['ERROR: "type" of sampling points was specified', ...
                    ' and arbitrary sample points are not supported!\n']); 
                return;
            end

            if ( isa(f_vals, 'function_handle') )
                %fprintf(2, 'Evaluating f_handle() at given points\n');
                [Obj.f_vals, Obj.op_dims] = Obj.i_pts.evaluate(f_vals); % Evaluate f at i_pts
            else
                % Check that the dimensions that have been passed in are correct
                s_fvals = size(f_vals);
                l_npts = length(Obj.i_pts.n_pts);
                if ( prod( s_fvals(end-l_npts+1:end) == Obj.i_pts.n_pts ) )
                    Obj.f_vals = f_vals;
                    Obj.op_dims = s_fvals(1:end-l_npts);
                else
                    error(['ERROR: Mismatch in dimensions of function values ', ...
                        'passed in and sample points\n']);
                end
            end

            Obj.in_dims = in_dims;

            % COLONS are automatically generated from the estimated dimension of the output
            Obj.colons = cell( size(Obj.op_dims) );
            [Obj.colons{:}] = deal(':');
            Obj.bounds = num2cell(bounds, 1);
            Obj.order = order;
        end

        function n_pts = getNPts(varargin)
            n_pts = varargin{1}.i_pts.getNPts(varargin{2:end});
        end

        function tpt = getPts(varargin)
            tpt = varargin{1}.i_pts.getPts(varargin{2:end});
        end

        function tpt = getPtAt(varargin)
            % This produces some answer (WRONG) for arbitrarily high positive integer values. FIX IT [TODO]
            tpt = varargin{1}.i_pts.getPtAt(varargin{2:end});
        end

        function der = firstDerivative(varargin)
            %{ DEBUG
            % fprintf(2, 'Received %d inputs. ', length(varargin));
            % fprintf(2, 'Computing derivative at: \n');
            % varargin
            % TODO: There was some extremely wierd bug here.
                % This doesn't work
                % [~, der] = varargin{1}.computeWithDer(varargin{2:end});

                % This, however, magically works
                % Obj = varargin{1};
                % args = {varargin{2:end}};
                % [~, der] = Obj.computeWithDer(args{:});
            %}
            Obj = varargin{1};
            args = {varargin{2:end}};
            [~, der] = Obj.computeWithDer(args{:});
        end

        function plotChebCoeffs(Obj)
            if ( prod(Obj.op_dims) > 3 )
                sub = num2cell(ind2sub(Obj.op_dims, 3)); % For a MOSFET device, 1 represents q(ds) and 3 represents i(ds)
            else
                sub = num2cell(ind2sub(Obj.op_dims, 1)); % For some use cases, like general functions, there may be just
                                                         % one output dimension and we might be looking for the plot for
                                                         % just that one
            end

            figure();
            if Obj.in_dims == 1
                series_coeffs = Obj.i_pts.getSeriesCoeffs( squeeze(Obj.f_vals(sub{:}, :))' );
                plot( abs(series_coeffs), 'LineWidth', 2.8 );
                xlabel('Coeff. Index', 'FontSize', 28);
                ylabel('Coeff. Value', 'FontSize', 28);
                set(gca, 'YScale', 'log');
            elseif Obj.in_dims == 2 
                series_coeffs = Obj.i_pts.getSeriesCoeffs( squeeze(Obj.f_vals(sub{:}, :, :))' );
                mesh( abs(series_coeffs) );
                set(gca, 'ZScale', 'log');
                xlabel('Coeff. Index_{i}', 'FontSize', 28);
                ylabel('Coeff. Index_{j}', 'FontSize', 28);
                ylabel('Coeff. Value(i,j)', 'FontSize', 28);
            else
                error('Too many dimensions for plotting series equivalent');
            end
            grid on;
            set(gca, 'FontSize', 28);
        end

        function plotDerivatives(Obj)
            ders = zeros(Obj.getNPts(1),1);
            for p_i = 1:length(ders)
                ders(p_i, 1) = Obj.firstDerivativeAtPt(p_i); 
            end
            figure, plot( Obj.getPts(1), ders, 'LineWidth', 2.0);
            grid on;
            set(gca, 'FontSize', 28);
        end

        function load(Obj, filename, prefix)
            Obj.i_pts = InterpolationPoints();
            Obj.i_pts.load([filename, '.sp']);
            % DEBUG
            % fprintf(2, 'Loading Interpolant object... \n');
            % END DEBUG
            load(filename, ...
                'in_dims', ...
                'op_dims', ...
                'order', ...
                'wts', ...
                'bounds', ...
                'colons', ...
                'extrap_slope', ...
                'f_vals', ...
                Obj.load_opts{:});
            Obj.in_dims         = in_dims;
            Obj.op_dims         = op_dims;
            Obj.order           = order;
            Obj.wts             = wts;
            Obj.bounds          = bounds;
            Obj.colons          = colons;
            Obj.extrap_slope    = extrap_slope;
            Obj.f_vals          = f_vals;
        end

        function save(Obj, filename, prefix)
            Obj.i_pts.save([filename, '.sp']);
            in_dims         = Obj.in_dims;
            op_dims         = Obj.op_dims;         
            order           = Obj.order;
            wts             = Obj.wts;             
            bounds          = Obj.bounds;
            colons          = Obj.colons;
            extrap_slope    = Obj.extrap_slope;
            f_vals          = Obj.f_vals;
            save(filename, ...
                'in_dims', ...
                'op_dims', ...
                'order', ...
                'wts', ...
                'bounds', ...
                'colons', ...
                'extrap_slope', ...
                'f_vals', ...
                Obj.save_opts{:});
        end
    end
end
