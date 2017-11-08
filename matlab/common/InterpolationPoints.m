classdef InterpolationPoints < SaveLoad
    properties (SetAccess = protected)
        order = 0;
        n_in_dims = 0;
        n_pts = [];
        bounds = [];    % Bounds in the global coordinates
        l_bounds = [];  % Bounds in the local coordinates
        pts = cell(0);
        sample_pts = cell(0);

        rs_factor = []; % Local variables required for change of coordinates
        center = [];
        n_center = [];
    end

    methods (Access = protected)
        function recenterAndSave(Obj, desired_pts)
            % This function gets a set of desired points and an
            % InterpolationPoints Obj. It then rescales the desired points
            % uniformly by a rescale factor that is decided by Obj.bounds and
            % uses the obtained sample_pts to sample the function that we
            % wish to interpolate

            % Reset the values of the various class properties
            Obj.n_pts =  zeros(1, Obj.n_in_dims);
            Obj.pts = cell(Obj.n_in_dims, 1);
            Obj.sample_pts = cell(Obj.n_in_dims, 1);
            Obj.l_bounds = Obj.bounds;
            Obj.rs_factor = ones(Obj.n_in_dims, 1);
            Obj.center = zeros(Obj.n_in_dims, 1);
            Obj.n_center = zeros(Obj.n_in_dims, 1);

            if ( nargin < 2 )   % Set everything to 0/1... 
                return;
            end

            if ( size(desired_pts, 1) ~= Obj.n_in_dims )
                error(['%d dimensional input data supplied. ', ...
                    '%d expected as a cell array. Aborting!\n'], ...
                    size(desired_pts,1), Obj.n_in_dims);
            end
        
            % Calculating class properties from the supplies points
            for d_i = 1:Obj.n_in_dims
                Obj.n_pts(1, d_i) = size(desired_pts{d_i}, 1);

                % First we squeeze these points into [-1 1]
                % TODO: While this works very well for Chebyshev polynomials,
                % which are defined in [-1, 1], it has some issues with DCT
                % based interpolant and should not be used there
                Obj.l_bounds(1, d_i) = desired_pts{d_i}(1);
                Obj.l_bounds(2, d_i) = desired_pts{d_i}(end);
                Obj.rs_factor(d_i, 1) = abs( Obj.bounds(2, d_i) - Obj.bounds(1, d_i) ) / ...
                    abs( Obj.l_bounds(2, d_i) - Obj.l_bounds(1, d_i) );
                Obj.center(d_i, 1) = ( Obj.l_bounds(2, d_i) + Obj.l_bounds(1, d_i) ) / 2.0;
                Obj.n_center(d_i, 1) = ( Obj.bounds(2, d_i) + Obj.bounds(1, d_i) ) / 2.0;
                Obj.sample_pts{d_i} = ( (desired_pts{d_i} - Obj.center(d_i, 1)) * ...
                    Obj.rs_factor(d_i, 1) ) + Obj.n_center(d_i, 1);
                Obj.pts{d_i} = desired_pts{d_i};
            end
        end
    end

    methods (Access = public)
        % Constructor for the class
        function obj = InterpolationPoints(n_in_dims, order, bounds, x_vals)
            if (nargin == 0)
                n_in_dims = 0;
                n_pts = [];
                %fprintf(2, 'WARNING: n_in_dims/n_pts not set!\n');
                return;
            end

            % Checking if the supplied n_pts matches the number of data dims
            if ( n_in_dims ~= size(order, 2) )
                error('Mismatch in n_in_dims and data dimensions.');
                return;
            end
            obj.n_in_dims = n_in_dims;
            obj.order = order;
            obj.bounds = zeros(2, n_in_dims);

            if (nargin > 2) % The bounds have been provided in an appropriately
                            % sized array

                % Check that we have values for each dimension
                bound_dim_check = size(bounds, 2);

                % Check that we have exactly 2 values (left, right) pair
                bound_lrp_check = size(bounds, 1);
                if (bound_dim_check ~= n_in_dims || bound_lrp_check ~= 2)
                    error(['Mismatch in bounds and data ', ...
                        'dimensions. Ensure bounds = 2 x n_in_dims'])
                    return;
                end

                obj.bounds = bounds;
                if (nargin > 3)
                    for d_i = 1:n_in_dims
                        d_pts =  x_vals{d_i};   % Has to be a cell array
                        if ( issorted(d_pts) )
                            obj.pts = x_vals;
                        else
                            error('Interpolation points not sorted');
                            return;
                        end
                    end
                end
            else    % Setting defaults [-1, 1] for all data dimensions
                obj.bounds(1,:) = -1;
                obj.bounds(2,:) = 1;
            end
        end 

        function [tpt, sub] = getPtAt(obj, pt, dim)
            if ( nargin > 2 )
                tpt = obj.sample_pts{dim}(pt);
                sub = {};
                return;
            end

            pt = int64(pt - 1);    % Special thanks to Matlab indexing
            tpt = zeros(obj.n_in_dims, 1);
            sub = cell(obj.n_in_dims, 1);  % Subscript index
            for d_i = 1:obj.n_in_dims
                sub{d_i} = 1 + rem(pt, obj.n_pts(d_i));
                tpt(d_i, 1) = obj.sample_pts{d_i}(sub{d_i});
                pt = idivide( pt, int64(obj.n_pts(d_i)) );
            end
        end

        function tpts = getPts(obj, dim)
            if (dim > obj.n_in_dims)
                error('ERROR: Index exceeds data dimensions');
                tpts = [];
                return; 
            end
            tpts = obj.sample_pts{dim};
        end

        function [f_vals, op_dims] = evaluate(Obj, f_handle)
            if ( ~strcmp(class(f_handle), 'function_handle') )
                error('Evaluation requested but supplied object ', ...
                    'is not a function handle');
                f_vals = [];
                return;
            end

            % Finding out the number of output dimensions for the function described by f_hanlde
            t_op = f_handle( zeros(Obj.n_in_dims, 1) );
            op_dims = size(t_op, 1);  % NOTE: Expect vector functions to return a column vector as an output.
                                            % Returning a row vector would cause it to be mistaken for a scalar function
                                            % and will definitiely lead to runtime errors

            f_vals = zeros([op_dims Obj.n_pts]);
            colons = cell(1, g_dims(op_dims));
            [colons{:}] = deal(':');

            for pt = 1:prod(Obj.n_pts)  % Can't think of a better way to do this
                                        % Even using ndgrid doesn't help much
                [tpt, sub] = Obj.getPtAt(pt);
                f_vals(colons{:}, sub{:}) = f_handle(tpt);
            end
        end

        function n_pts = getNPts(Obj, d_i)
            % This function returns the number of sample points that have been embedded in the object. If a dimension
            % index is supplied, the number of sample points along that dimension are returned, otherwise, the total
            % number of sample points is returned
            if ( nargin > 1 )
                n_pts = Obj.n_pts(d_i);
            else
                n_pts = prod(Obj.n_pts);
            end
        end

        function pt_index = findPt(Obj, pt)
            % Get the indices corresponding to the supplied "pt" in the stored cell array of points. Keep in mind that
            % Chebyshev points are stored in a descending order.
            pt_index = zeros(1, Obj.n_in_dims);
            for d_i = 1:Obj.n_in_dims
                p_diff = (pt >= Obj.getPts(d_i));
                % Since Obj.pts{ ... } are arranged in a descending order, pt - Obj.pts{ ... } should be arranged in an
                % ascending order. This lets us use MATLAB's builtin function to find it
               pt_index(1, d_i) = find(p_diff, 1); 
            end
        end

        function [x_out, dx_out] = rescaleShiftInput(Obj, x_in)
            x_out = ( (x_in - Obj.n_center) ./ Obj.rs_factor ) + Obj.center;
            dx_out = 1 ./ Obj.rs_factor;
        end

        function load(Obj, filename)
            % DEBUG
            % fprintf(2, 'Loading Interpolantion Points...\n');
            % END DEBUG
            load(filename,  ...
                'order', ...
                'n_in_dims', ...
                'n_pts', ...
                'bounds', ...
                'l_bounds', ...
                'pts', ...
                'sample_pts', ...
                'rs_factor', ...
                'center', ...
                'n_center', ...
                Obj.load_opts{:});
            Obj.order       = order;
            Obj.n_in_dims   = n_in_dims;
            Obj.n_pts       = n_pts;
            Obj.bounds      = bounds;
            Obj.l_bounds    = l_bounds;
            Obj.pts         = pts;
            Obj.sample_pts  = sample_pts;
            Obj.rs_factor   = rs_factor;
            Obj.center      = center;
            Obj.n_center    = n_center;
        end

        function save(Obj, filename)
            order       = Obj.order;
            n_in_dims   = Obj.n_in_dims;
            n_pts       = Obj.n_pts;
            bounds      = Obj.bounds;
            l_bounds    = Obj.l_bounds;
            pts         = Obj.pts;
            sample_pts  = Obj.sample_pts;
            rs_factor   = Obj.rs_factor;
            center      = Obj.center;
            n_center    = Obj.n_center;
            save(filename,  ...
                'order', ...
                'n_in_dims', ...
                'n_pts', ...
                'bounds', ...
                'l_bounds', ...
                'pts', ...
                'sample_pts', ...
                'rs_factor', ...
                'center', ...
                'n_center', ...
                Obj.save_opts{:});
        end
    end
end
