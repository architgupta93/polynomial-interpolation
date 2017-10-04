% Class definitions for an interpolant object
classdef Interpolant < SaveLoad
    properties (SetAccess = protected)
        n_in_dims = 0;
        op_dims = 0;
        order = [];
        wts = [];
        i_pts = InterpolationPoints();
        bounds = cell(0);
        colons = cell(0);
        extrap_slope = [];
        f_vals = [];
    end

    methods (Static = true, Access = protected)
        function success = isColumnVector(x)
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
            var = [];
            der = [];
        end
        function der = secondDerivative(Obj, x_in);
            val = [];
            der = [];
        end
    end

    methods (Access = public)
        function Obj = Interpolant(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)
            if (nargin == 0)
                % DEBUG MESSAGE [TODO]: Remove this later on
                % fprintf(2, 'Instantiating an EMPTY interpolant of class: %s! \n', class(Obj)); 
                return;
            end

            % Checking if interpolation type (i_type) has been provided
            if ( isstr(i_type_or_x_vals) )
                if ( strcmp(i_type_or_x_vals, 'uniform') || ...
                    strcmp(i_type_or_x_vals, 'u') )
                    Obj.i_pts = UniformPoints(n_in_dims, order, bounds);
                elseif ( strcmp(i_type_or_x_vals, 'chebyshev') || ...
                    strcmp(i_type_or_x_vals, 'c') )
                    Obj.i_pts = ChebyshevPoints1Analyst(n_in_dims, order, bounds);
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

            if ( strcmp( class(f_vals), 'function_handle' ) )
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

            Obj.n_in_dims = n_in_dims;
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
            if Obj.n_in_dims == 1
                series_coeffs = Obj.i_pts.getSeriesCoeffs( squeeze(Obj.f_vals(sub{:}, :))' );
                plot( abs(series_coeffs), 'LineWidth', 2.8 );
                xlabel('Coeff. Index', 'FontSize', 28);
                ylabel('Coeff. Value', 'FontSize', 28);
                set(gca, 'YScale', 'log');
            elseif Obj.n_in_dims == 2 
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
                'n_in_dims', ...
                'op_dims', ...
                'order', ...
                'wts', ...
                'bounds', ...
                'colons', ...
                'extrap_slope', ...
                'f_vals', ...
                Obj.load_opts{:});
            Obj.n_in_dims       = n_in_dims;
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
            n_in_dims       = Obj.n_in_dims;
            op_dims         = Obj.op_dims;         
            order           = Obj.order;
            wts             = Obj.wts;             
            bounds          = Obj.bounds;
            colons          = Obj.colons;
            extrap_slope    = Obj.extrap_slope;
            f_vals          = Obj.f_vals;
            save(filename, ...
                'n_in_dims', ...
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
