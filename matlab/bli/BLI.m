classdef BLI < Interpolant
    properties (SetAccess = protected)
        % Store the coefficients for the numerator and denominator separately
        nu_vals = []; % Coefficients for the numerator
        de_vals = []; % Coefficients for the denominator
    end

    methods (Access = public)
        function Obj = BLI(varargin)
        % function Obj = BLI(f_vals, in_dims, bounds, order, i_type_or_x_vals);
        % Class constructor
            Obj = Obj@Interpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            if ( iscell(varargin{end}) )
                Obj.coeffs = tensorProduct(Obj.coeffs);
            else
                Obj.coeffs = tensorProduct(Obj.i_pts.getBLIWeights());
            end
            diff_dims = ndims(Obj.coeffs) - ndims(Obj.f_vals); 
            Obj.nu_vals = times(Obj.f_vals, shiftdim(Obj.coeffs', diff_dims));
            Obj.de_vals = Obj.coeffs';

            Obj.extrap_slope = zeros([Obj.op_dims 2]);
            Obj.extrap_slope(Obj.colons{:}, 1) = Obj.firstDerivativeAtPt(1);
            Obj.extrap_slope(Obj.colons{:}, 2) = Obj.firstDerivativeAtPt(); % Calculates the value at the end
        end

        function der = firstDerivativeAtPt(Obj, pt_index)
            if ( nargin < 2)
                pt_index = Obj.i_pts.n_pts(1);
            end

            all_pts = Obj.i_pts.pts{1};
            this_pt = all_pts(pt_index);
            exclude_pts = 1 ./ (this_pt - all_pts);
            exclude_pts(pt_index) = 0;  % Removing this entry from the multiplication... NOTE HERE:
                                        % Quite sure that all the unnecessary copying would lead to more resource
                                        % consumption than multiplication of 1 row by 0 and subsequent addition
            diff_dims = ndims(exclude_pts) - ndims(Obj.nu_vals); 

            der = ( -times( Obj.f_vals(:, pt_index), Obj.de_vals * exclude_pts ) + ...
               sum(times(Obj.nu_vals, shiftdim(exclude_pts', diff_dims)), ndims(Obj.nu_vals)) ) / Obj.de_vals(pt_index); 

            % Previous implementation: TODO If anything goes wrong, uncomment this
            %exclude_pts = 1 ./ (this_pt - [ all_pts(1:pt_index-1, 1); all_pts(pt_index+1:end, 1) ]);
            % exclude_nu = [Obj.nu_vals(:, 1:pt_index-1) Obj.nu_vals(:, pt_index+1:end)];
            % exclude_de = [Obj.de_vals(1, 1:pt_index-1) Obj.de_vals(1, pt_index+1:end)];
            %der = ( - ( Obj.f_vals(:, pt_index) * exclude_de * exclude_pts ) + ...
            %   (exclude_nu * exclude_pts) ) / Obj.de_vals(pt_index); 
        end

        function der = secondDerivativeAtPt(Obj, pt_index)
            der = 0;
            %TODO: Will be required for better smoothing
        end

        function [val, der] = computeWithDer(Obj, x_in, coeffs)
            %{ DEBUG
            %    Obj
            %    x_in
            %}
            [x_in, dx_out] = Obj.i_pts.rescaleShiftInput(x_in);
            if (nargin < 3)
                coeffs = Obj.nu_vals;
            end

            x_vals = Obj.i_pts.pts{1}; 
            x_rel = x_vals' - x_in;  % x_in is a scalar and x_vals is a vector

            % Some pre-emptive code for extrapolation... Maybe we can do this more neatly [TODO]
            % [TODO]: This needs to be done properly for the case where coefficients are supplied externally. In that
            % case, the extrapolation slope is no longer valid
            if ( x_rel(1) <= 0.0 )   % Remeber getPts returns values in decreasing order
                der = dx_out * Obj.extrap_slope(Obj.colons{:}, 1);
                val = rdivide(coeffs(Obj.colons{:}, 1), Obj.de_vals(1)) - x_rel(1) * der;
                return;
            elseif ( x_rel(end) >= 0.0 )
                der = dx_out * Obj.extrap_slope(Obj.colons{:}, 2);
                val = rdivide(coeffs(Obj.colons{:}, length(x_rel)), Obj.de_vals(length(x_rel))) - x_rel(end) * der;
                return;
            end

            % If x_in is exactly equal to one of the points in i_pts, return the
            % corresponding value. We already know that the values in x_rel are
            % sorted. So, life is good. Use matlab's MEX function: Careful or
            % the function will crash because of segmentation fault

            %[found, loc] = builtin('_ismemberhelper', 0, x_rel);
            %[HACK] Just to get the timing down
            found = 0; loc = 1;
            if (found)
                % Too many infinities to deal with, we will just use the formula that we already have
                val = rdivide(coeffs(Obj.colons{:}, loc), Obj.de_vals(loc));
                der = Obj.fitrstDerivativeAt(loc);
            else
                diff_dims = ndims(x_rel) - ndims(coeffs); 
                l_nu_val = sum ( rdivide(coeffs, reshape(x_rel, [ones(1, abs(diff_dims)), size(x_rel)])), ndims(coeffs));
                l_de_val = sum ( rdivide(Obj.de_vals, x_rel), 2);
                val = rdivide(l_nu_val, l_de_val);

                if (nargout > 1)
                    % Welcome to a new world of pain: It is too difficult to even try to write the whole expression here.
                    % Please refer to notes "Working-out-continuity-for-piecewise-Chebyshev-series" for details. The
                    % terminology here has also been borrowed from the expressions there, so it will be hard to follow
                    % anything if they are not readily available

                    % Calculating j_sum (for each i first)
                    x_i_minus_j = x_rel' - x_rel; 
                    
                    % Keep in mind that x_rel is a row vector. This means that x_rel - x_rel' = xj - xi,  whereas, doing
                    % x_rel' - x_rel = xi - xj. Since we will be pre-multiplying with a row-vector, we want the latter to
                    % happen, i.e., the quantity that is varying across the rows should carry a positive sign.
                    % Therefor, the substraction means that x_i_minus_j has i along the columns and j along the rows 

                    j_vals = (Obj.de_vals ./ ( x_rel .^ 2 ));
                    j_sum_times_xi_sq = ( j_vals * x_i_minus_j) ./ ( x_rel .^ 2 );  

                    % Obj.de_vals is a row vector. 
                    % Obj.de_vals ./ ( x_rel .^ 2 ) corresponds to wj/(x-x_i)^2 term which appears both in the numerator and
                    % the denominator in the overall expression. In the numerator, these are weighted by (x_i - x_j)

                    l_nu_der = sum ( times(coeffs, reshape(j_sum_times_xi_sq, [ones(1, abs(diff_dims)), size(x_rel)])), ndims(coeffs));
                    l_de_der = l_de_val .^ 2;
                    der = dx_out * rdivide(l_nu_der, l_de_der);
                end
            end
        end

        function der = secondDerivative(obj, x_in)
            der = 0;
        end

        function load(Obj, filename)
            % Call a load on the superclass
            load@Interpolant(Obj, filename);

            % Load the remaining class members in an add-on file
            % DEBUG
            % fprintf(2, 'Loading BLI Add-ons...\n');
            % END DUBUG
            load([filename, '.bliao'], ... 
                'nu_vals', ...
                'de_vals', ...
                Obj.load_opts{:});
            Obj.nu_vals = nu_vals;
            Obj.de_vals = de_vals;
        end

        function save(Obj, filename)
            % Call a save on the superclass
            save@Interpolant(Obj, filename);

            % Save the remaining class members as an add-on
            % DEBUG
            % fprintf(2, 'Saving BLI Add-ons...\n');
            % END DEBUG
            nu_vals = Obj.nu_vals;
            de_vals = Obj.de_vals;
            save([filename, '.bliao'], ...
                'nu_vals', ...
                'de_vals', ...
                Obj.save_opts{:});
        end
    end
end
