classdef PiecewiseInterpolant < SaveLoad
    properties (Constant, Access = protected)
        d_order = 6;
    end

    properties (SetAccess = protected, GetAccess = public)
        n_in_dims = 0;     % Number of input dimensions
        op_dims = 0;
        n_pieces = [];  % Array containing numer of pieces in each dimension 
        order = 0;
        is_smooth = false;
        bounds = cell(0);
        colons = cell(0);   % Somewhat useful for accesses
    end

    properties (Access = protected)
        m_interp = cell(0);
    end

    methods (Access = protected) % Should actually be abstract functions but that isn't supported very well in Octave
        function populate(Obj)
            Obj.m_interp = squeeze(cell([Obj.n_pieces 1]));
            for interp_in = 1 : prod(Obj.n_pieces)
                Obj.m_interp{interp_in} = Interpolant();
            end
            if ( isempty(Obj.n_pieces) )
                warning('Running populate with Abstract class "PiecewiseInterpolant" while n_pieces is empty. Hopefully you know what is going on :)');
            end
        end
    end

    methods (Access = protected)
        function l_bounds = getLocalBounds(Obj, pc_index)
            b_i = pc_index - 1;
            l_bounds = zeros(2, Obj.n_in_dims);
            for d_i = 1:Obj.n_in_dims
                l_index = rem(b_i, Obj.n_pieces(d_i));
                l_bounds(1, d_i) = Obj.bounds{d_i}(l_index+1); 
                l_bounds(2, d_i) = Obj.bounds{d_i}(l_index+2);
                b_i = idivide( int64(b_i), int64(Obj.n_pieces(d_i)) );
            end
        end

        function pc = getPieceIndex(Obj, x_in)
            pc = cell(Obj.n_in_dims, 1);
            for d_i = 1:Obj.n_in_dims
                pc{d_i, 1} = find(x_in(d_i, 1) <= [ -Inf; Obj.bounds{d_i}(2:end-1); Inf ], 1) - 1;

                % For debugging purposes only
                if (pc{d_i, 1} < 1)
                    error('Reported index of input %d is %d', x_in(d_i), pc{d_i, 1});
                end
            end
        end

        function o_qty = matchDimsWithNPieces(Obj, qty)
            % Checking that order is appropriately sized, i.e., either a
            % scalar value denoting the required number of points in each
            % dimension OR the same size as n_pieces, i.e., telling how many
            % points to use in each piece.
            if ( length(qty) == 1 )
                o_qty = qty * ones( size(Obj.n_pieces) );
                fprintf(2, 'WARNING: Supplied SCALAR qty. Expanding to a vector.\n');
            elseif ( prod(size(qty) == size(Obj.n_pieces)) )
                o_qty = qty;
            else
                error(['ERROR: Supplied vector for quantity', ...
                    'inappropriately sized w.r.t. n_pieces\n']);
            end
        end

        function [val, der] = interpolate(Obj, x_in)
            pc_index = Obj.getPieceIndex(x_in);
            if (isempty([pc_index{:}]))
                error('Piecewise Interpolant: Failed to provide queried interpolant value @ [%d] ', full(x_in));
            end

            % Saving some computation time if the derivative is not required
            if (nargout > 1)
                [val, der] = Obj.m_interp{pc_index{:}}.computeWithDer(x_in);
            else
                val = Obj.m_interp{pc_index{:}}.computeWithDer(x_in);
            end
        end
    end

    methods (Access = public)  % Abstract functions from MATLAB
    % API functions that all the inheriting classes should implement
        function [val, der] = computeWithDer(Obj, x_in);
            val = [];
            der = [];
        end
    end

    methods (Access = public)
        function Obj = PiecewiseInterpolant(f_handle, n_in_dims, bounds, order, i_type, smooth)
            if ( nargin == 0 )
                % DEBUG
                fprintf(2, 'Instantiating an EMPTY Piecewise Interpolant\n');
                return;
            end
            Obj.n_in_dims = n_in_dims;

            % Either we are given a function handle and we can use it to  determine the dimensions of the output that
            % the function produces OR we could be given the function values (appropriately sized) and we would have to
            % figure out the rest from there
            if ( strcmp(class(f_handle), 'function_handle') )
                tval = f_handle(zeros(n_in_dims, 1));
                Obj.op_dims = g_size(tval);
            elseif ( isnumeric(f_handle) )
                % We have some GIANT matrix that needs to be resolved. However, this matrix is arranged as follows:
                % data x n_pieces x n_pts (per piece)
                s_fvals = size(f_handle);
                Obj.op_dims = s_fvals(1:end - ( Obj.n_in_dims +   1 ));
                %                                   ^^              ^^
                %                               Accounting for   n_pieces
                %                                the input(s)       (1)
            else
                error('ERROR: Invalid/Unsupported data type for input function\n');
            end
            Obj.colons = cell(size(Obj.op_dims));
            [Obj.colons{:}] = deal(':');

            if ( isa(bounds, 'cell') )
                Obj.bounds = bounds;
            else
                error('Expecting a cell array for bounds, indexing each dimension in the columns');
            end

            Obj.n_pieces = zeros( size(Obj.bounds) );
            for pc_index = 1 : size(Obj.n_pieces, 2)
                Obj.n_pieces(1, pc_index) = size( Obj.bounds{pc_index}, 1 ) - 1;
            end

            if ( size(Obj.n_pieces,2) ~= n_in_dims )
                fprintf(2, 'NOTE: Expect n_pieces to be a row vector\n');
                error(['ERROR: Mismatch in supplied n_pieces ', ...
                    '(per dimensions) and n_in_dims\n']);
            end

            if (nargin > 2)
                % Checking that the supplied bounds match the intended number of
                % pieces for interpolation using piecewise Chebyshev series
                if ( (size(bounds, 2) ~= Obj.n_in_dims) )
                    error(['ERROR: Mismatch in supplied bounds and ', ...
                        'n_in_dims. Expecting 1 x n_in_dims cell array\n']);
                    return;
                else
                    for d_i = 1:Obj.n_in_dims
                        if ( size(bounds{d_i}) ~= (Obj.n_pieces(d_i) + 1) )
                            fprintf(2, ['ERROR: Mismatch in supplied bounds and ', ...
                                'n_pieces at dimension %s!\n'], d_i);
                            return;
                        end
                    end
                end
            else
                % We will take a n_in_dims dimensional cube (x[-1, 1])^n_in_dims and
                % divide it into n_pieces equal pieces "TODO"
                fprintf(2, ['ERROR: Automatic division not supported ', ... 
                    'at the moment\n']);
            end

            if (nargin > 3)
                Obj.order = Obj.matchDimsWithNPieces(order);
            else
                % Use the default value for order for each piece (each
                % dimensiona as well) supplied in the properties (Constant)
                % section
                Obj.order = Obj.d_order * ones( size(Obj.n_pieces) );
            end

            if (nargin > 5)
                Obj.is_smooth = true;
                % Otherwise, use the default value, which is FALSE
            end
        end

        function der = secondDerivative(Obj, x_in)
            pc_index = Obj.getPieceIndex(x_in);
            der = Obj.m_interp{pc_index{:}}.secondDerivative(x_in);
        end

        function der = firstDerivative(Obj, x_in)
            [~, der] = Obj.computeWithDer(x_in);
        end

        function plotChebCoeffs(Obj, pc_index)
            Obj.m_interp{pc_index}.plotChebCoeffs();
        end

        function plot(Obj, pc_index)
            if ( nargin < 2 )
                pc_index = cell(Obj.n_in_dims,1);
                for d_i = 1:Obj.n_in_dims
                    pc_index{d_i, 1} = 1;
                end
            end

            if ( prod(pc_index{:}) > prod( size(Obj.m_interp) ) || ...
                 prod(pc_index{:}) < 1 )
                fprintf(2, ['ERROR: No interpolant pieces found for ', ...
                    'plotting OR index exceeds the data dimensions\n']);
            else
                Obj.m_interp{pc_index{:}}.plot();
            end
        end

        function save(Obj, dirname)
            % Save a Piecewise interpolant object. This function creates a new directory and puts all the associated
            % files in that directory.
            % [TODO]: Check if the directory already exists and throw a warning asking the user if he/she wants to
            % overwrite the existing save and create a new directory in its place with the same name.
            if ( exist(dirname, 'dir') )
                waiting_for_yes_no = true; 
                while ( waiting_for_yes_no )
                    yes_to_overwrite = input(['Data directory ', dirname, ' exists! Want to overwrite? [y/n]'], 's');
                    if (yes_to_overwrite == 'y' || yes_to_overwrite == 'Y')
                        rmdir(dirname,'s');
                        waiting_for_yes_no = false;
                    elseif (yes_to_overwrite == 'n' || yes_to_overwrite == 'N')
                        fprintf(2, 'Aborting save on user request\n');
                        return;
                    end
                end
            end

            mkdir(dirname);
            % fprintf(2, 'Saving Piecewise interpolant...\n');
            % fprintf(2, 'Saving Header Information...\n');
            % Saving class members that aren't instances of other classes themselves
            n_in_dims   = Obj.n_in_dims;
            op_dims     = Obj.op_dims;
            n_pieces    = Obj.n_pieces;
            order       = Obj.order;
            is_smooth   = Obj.is_smooth;
            bounds      = Obj.bounds;
            colons      = Obj.colons;
            save([dirname, '/header.info'], ...
                'n_in_dims', ... 
                'op_dims', ...
                'n_pieces', ...
                'order', ...
                'is_smooth', ...
                'bounds', ...
                'colons', ...
                Obj.save_opts{:});   

            % Using linear indexing for picking the individual interpolants for saving them
            % fprintf(2, 'Saving Interpolants...\n');
            for interp_in = 1 : prod(Obj.n_pieces)
                Obj.m_interp{interp_in}.save([dirname, '/id__', num2str(interp_in), '.intrp']);
            end
            % fprintf(2, '\n');
        end

        function load(Obj, dirname)
            % Loading class members that aren't instances of other classes themselves. Refer to the save function if it
            % to see which file are these located in. Current location is filename
            % fprintf(2, 'Loading Piecewise interpolant...\n');
            % fprintf(2, 'Loading Header Information...\n');
            load([dirname, '/header.info'], ...
                'n_in_dims', ... 
                'op_dims', ...
                'n_pieces', ...
                'order', ...
                'is_smooth', ...
                'bounds', ...
                'colons', ...
                Obj.load_opts{:});   
            Obj.n_in_dims   = n_in_dims;
            Obj.op_dims     = op_dims;
            Obj.n_pieces    = n_pieces;
            Obj.order       = order;
            Obj.is_smooth   = is_smooth;
            Obj.bounds      = bounds;
            Obj.colons      = colons;

            % Each child class has to implement a populate function. Which now fills the m_interp array with instance of
            % the appropriate class
            Obj.populate();

            % Again using linear indexing
            % fprintf(2, 'Loading Interpolants...\n');
            for interp_in = 1 : prod(Obj.n_pieces)
                interp_filename = [dirname, '/id__', num2str(interp_in), '.intrp'];
                % Checking if the file exists
                if ( ~exist(interp_filename) )
                    error( 'Interpolant ID: %d expected but NOT FOUND in file %s', interp_in, interp_filename );
                end
                Obj.m_interp{interp_in}.load(interp_filename);
                % fprintf(2, 'Loaded %s!\n\n', interp_filename);
            end
        end
    end
end