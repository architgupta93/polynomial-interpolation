classdef PiecewiseInterpolant < SaveLoad
    properties (Constant, Access = protected)
        d_order = 6;
    end

    properties (SetAccess = protected, GetAccess = public)
        % INPUT: As of now, we only allow functions of vectors to be
        % interpolated, a scalar number, therefore, defines its size
        in_dims = 0;

        % OUTPUT: This can have an arbitrary dimension, as a result, the
        % dimensionality needs to be stored as a row vector
        op_dims   = 0;

        % N PIECES: The number of pieces per input dimension. Only a Cartesian
        % partitioning of the overall input space is allowed with piecewise
        % inteprolants.
        n_pieces  = [];

        % REFINENESS: For each interpolant, there is a notion of 'refinement'.
        % ORDER defines the degree of refineness. Higher the order, larger the
        % size of the interpolant, and hopefully, better the interpolation. The
        % exact nature of this depends on the interpolant
        % TODO: Presently, all the pieces in a dimension are assigned the same
        % order. This is not necessary!
        order     = 0;

        % IS SMOOTH: Do we smoothen the pieces using an overlaid interpolant?
        is_smooth = false;

        % DOMAIN: Domain bounds for each piece inside which the interpolant CAN
        % interpolate the values. Each cell has to define the boundary values
        % for a dimension.
        bounds    = cell(0);

        % I_TYPE: The type of interpolant (sample points) to be used
        i_type    = [];

        % COLONS: This is a MATLAB/OCTAVE centric field required for indexing
        % all the OUTPUTs corresponding to a specific (partial) input value.
        % For example, if the output has dimesions 2x3x4x5, then for any input
        % value, these values are indexed using {:, :, :, :}, which is stored
        % in colons
        colons    = cell(0);

        % ACCESS_HANDLE: A Function handle which has to be a constructor for the
        % interpolant class. This can be used to populate the piecewise
        % interpolant
        acc_han   = @Interpolant;
    end

    properties (Access = protected)
        m_interp = cell(0);
    end

    methods (Access = protected)
        function populate(Obj, f_handle)
            Obj.m_interp = squeeze(cell([Obj.n_pieces 1]));
            if (nargin > 1)
                for interp_in = 1 : prod(Obj.n_pieces)
                    % This should now have all the required fields
                    l_bounds                = {Obj.getLocalBounds(interp_in)};
                    Obj.m_interp{interp_in} = Obj.acc_han(f_handle, Obj.in_dims, ...
                        l_bounds, Obj.order, Obj.i_type);
                    % After the interpolant has been fit, smoothing needs to be
                    % done by supplying the derivative and value for the
                    % adjoining interpolants
                    % TODO: This needs to be changed for a generalized spline
                    if (Obj.is_smooth && interp_in > 1)
                        corner_slope = Obj.m_interp{interp_in}.firstDerivativeAtPt();
                        fprintf('DEBUG: Fitting corner slope: %d\n', ...
                            corner_slope);
                        Obj.m_interp{interp_in-1}.fit([], corner_slope);
                    end
                end
            end

            if ( isempty(Obj.n_pieces) )
                warning('Running populate with Abstract class "PiecewiseInterpolant" while n_pieces is empty. Hopefully you know what is going on :)');
            end
        end
    end

    methods (Access = protected)
        function l_bounds = getLocalBounds(Obj, pc_index)
            b_i = pc_index - 1;
            l_bounds = zeros(2, Obj.in_dims);
            for d_i = 1:Obj.in_dims
                l_index = rem(b_i, Obj.n_pieces(d_i));
                l_bounds(1, d_i) = Obj.bounds{d_i}(l_index+1); 
                l_bounds(2, d_i) = Obj.bounds{d_i}(l_index+2);
                b_i = idivide( int64(b_i), int64(Obj.n_pieces(d_i)) );
            end
        end

        function pc = getPieceIndex(Obj, x_in)
            pc = cell(Obj.in_dims, 1);
            for d_i = 1:Obj.in_dims
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
                error(['ERROR: Supplied vector for quantity inappropriately sized w.r.t. n_pieces.', ...
                    ' Expecting same size as n_pieces']);
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

        function setAccessHandle(Obj, access_handle, f_handle)
            Obj.acc_han = access_handle; 
            Obj.populate(f_handle);
        end
    end

    methods (Access = public)  % Abstract functions from MATLAB
    % API functions that all the inheriting classes should implement
        function [val, der] = computeWithDer(Obj, x_in);
        % function [val, der] = computeWithDer(Obj, x_in);
        % DEFAULT IMPLEMENTATION
            if (nargout > 1)
                [val, der] = Obj.interpolate(x_in);
            else
                val = Obj.interpolate(x_in);
            end
        end
    end

    methods (Access = public)
        function Obj = PiecewiseInterpolant(varargin)
        % function Obj = PiecewiseInterpolant(f_handle, in_dims, bounds, ...
        %    order, i_type, access_handle, smooth)
        % Class constructor
        % Takes in the following input(s):
        %   F_HANDLE: Function handle (or values) for the function to be
        %       interpolated.
        %   IN_DIMS: Input dimensions (Number).
        %   BOUNDS: CELL ARRAY of the breakpoints in each dimension.
        %   ORDER: See REFINENESS in the description of properties above.
        %   I_TYPE: The category of sample points to be used.
        %   ACCESS_HANDLE: A Class constructor handle to be used to fill in the
        %       individual polynomial pieces
        %   SMOOTH: [BOOL] Whether the individual polynomial interpolants
        %       should be smoothened using and overlaid interpolant.

            if ( nargin < 2 )
                % DEBUG
                fprintf(2, 'Instantiating an EMPTY Piecewise Interpolant\n');
                return;
            end
            f_handle    = varargin{1};
            Obj.in_dims = varargin{2};

            % Either we are given a function handle and we can use it to
            % determine the dimensions of the output that the function produces
            % OR we could be given the function values (appropriately sized) and
            % we would have to figure out the rest from there
            if ( strcmp(class(f_handle), 'function_handle') )
                % Get a test value, use it to determine the output size
                tval = f_handle(zeros(in_dims, 1));
                Obj.op_dims = g_size(tval);
            elseif ( isnumeric(f_handle) )
                % We have some GIANT matrix that needs to be resolved. However,
                % this matrix is arranged as follows: data x n_pieces x n_pts
                % (per piece)
                % TODO: Something has to be done here! Not sure whaT
                s_fvals = size(f_handle);
                Obj.op_dims = s_fvals(1:end - ( Obj.in_dims +   1 ));
                %                                   ^^              ^^
                %                               Accounting for   n_pieces
                %                                the input(s)       (1)
            else
                error('Invalid/Unsupported data type for input function');
            end
            Obj.colons = cell(size(Obj.op_dims));
            [Obj.colons{:}] = deal(':');

            if ( isa(bounds, 'cell') )
                Obj.bounds = varargin{3};
            else
                error('Expecting a cell array for bounds, indexing each dimension in the columns');
            end

            Obj.n_pieces = zeros( size(Obj.bounds) );
            for pc_index = 1 : size(Obj.n_pieces, 2)
                Obj.n_pieces(1, pc_index) = size( Obj.bounds{pc_index}, 1 ) - 1;
            end

            if ( size(Obj.n_pieces,2) ~= in_dims )
                fprintf(2, 'NOTE: Expect n_pieces to be a row vector\n');
                error('Mismatch in supplied n_pieces (per dimensions) and in_dims');
            end

            if (nargin > 2)
                % Checking that the supplied bounds match the intended number of
                % pieces for interpolation using piecewise Chebyshev series
                if ( (size(bounds, 2) ~= Obj.in_dims) )
                    error(['Mismatch in supplied bounds and ', ...
                        'in_dims. Expecting 1 x in_dims cell array\n']);
                    return;
                else
                    for d_i = 1:Obj.in_dims
                        if ( size(bounds{d_i}) ~= (Obj.n_pieces(d_i) + 1) )
                            error(['Mismatch in supplied bounds and ', ...
                                'n_pieces at dimension %s!'], d_i);
                            return;
                        end
                    end
                end
            else
                % We will take a in_dims dimensional cube (x[-1, 1])^in_dims and
                % divide it into n_pieces equal pieces "TODO"
                error(['Automatic division not supported at the moment']);
            end

            if (nargin > 3)
                Obj.order = Obj.matchDimsWithNPieces(varargin{4});
            else
                % Use the default value for order for each piece (each
                % dimensiona as well) supplied in the properties (Constant)
                % section
                Obj.order = Obj.d_order * ones( size(Obj.n_pieces) );
                return;
            end

            if (nargin > 4)
                Obj.i_type = varargin{5};
            else
                Obj.i_type = 'uniform';
                return;
            end

            if (nargin > 5) 
                % We have the access handle to the class constructor which is an
                % interpolant. We can use this to construct the piecewise
                % interpolant
                Obj.acc_han = varargin{6};
                Obj.populate(f_handle);
            else
                return;
            end

            if (nargin > 6)
                Obj.is_smooth = varargin{7};
                % Otherwise, use the default value, which is FALSE
            end
        end

        function clearAllInterpolants(Obj)
            % Clear the contents of m_interp so that a new interpolant can be
            % put in
            Obj.m_interp = {};
        end

        function setInterpolant(Obj, pc_idx, access_handle, f_handle, ...
            bounds, order, i_type)
        % function setInterpolant(Obj, pc_idx, access_handle, f_handle, ...
        %     bounds, order, i_type)
        % Replace an existing interpolant with one of your choice. Optionally,
        % you can also pick the parameters for this new interpolant
            if (nargin < 7)
                i_type = Obj.i_type;
                if (nargin < 6)
                    order = Obj.order;
                    if (nargin < 5)
                        bounds = {Obj.getLocalBounds(pc_idx)};
                    end
                end
            end

            Obj.m_interp{pc_idx} = access_handle(f_handle, Obj.in_dims, bounds, ...
                order, i_type);
        end

        function ironOut(Obj, pc_idx)
            r_der = [];
            l_der = [];

            % See if right derivative needs to be matched
            if (pc_idx < length(Obj.m_interp))
                r_der = Obj.m_interp{pc_idx+1}.firstDerivativeAtPt();
            end

            % See if left derivative needs to be matched
            if (pc_idx > 1)
                l_der = Obj.m_interp{pc_idx-1}.firstDerivativeAtPt(1);
            end

            % Match the specified derivatives
            Obj.m_interp{pc_idx}.fit(l_der, r_der);
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
                pc_index = cell(Obj.in_dims,1);
                for d_i = 1:Obj.in_dims
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
            in_dims   = Obj.in_dims;
            op_dims     = Obj.op_dims;
            n_pieces    = Obj.n_pieces;
            order       = Obj.order;
            is_smooth   = Obj.is_smooth;
            bounds      = Obj.bounds;
            colons      = Obj.colons;
            save([dirname, '/header.info'], ...
                'in_dims', ... 
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

        % TODO: Save/Load functions need to be updated after all the changes
        % that have been made to the class.
        function load(Obj, dirname)
            % Loading class members that aren't instances of other classes themselves. Refer to the save function if it
            % to see which file are these located in. Current location is filename
            % fprintf(2, 'Loading Piecewise interpolant...\n');
            % fprintf(2, 'Loading Header Information...\n');
            load([dirname, '/header.info'], ...
                'in_dims', ... 
                'op_dims', ...
                'n_pieces', ...
                'order', ...
                'is_smooth', ...
                'bounds', ...
                'colons', ...
                Obj.load_opts{:});   
            Obj.in_dims   = in_dims;
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
