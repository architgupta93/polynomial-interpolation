classdef Spline1D < SplineInterpolant
    methods (Access = public)
        function Obj = Spline1D(varargin)
        % function Obj = Spline1D(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)
        % Class Constructor
            Obj = Obj@SplineInterpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            xi = Obj.i_pts.pts{1};
            n_plus_1 = Obj.i_pts.getNPts(1);

            if ( n_plus_1 <= 3)
                error(['You can write your interpolation on a piece of paper for ', ...
                    'so few points. Go get a life!']);
            else
                % Now that we have more than 3 points to work with, we will set up
                % sparse matrix equation that gives us the coefficients of the spline
                % interpolation 

                h = diff(xi)';   % Get the deltas between the individual xi elements
                % DIFF(X), for a vector X, is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)].
                if (~( isempty(find(h<=0,1)) || isempty(find(h>=0,1)) ))
                %% Checking that h DOES NOT contain both negative and positive entries
                %% the entries must be sorted (In increasing or decreasing order)
                    error('ERROR: Xi should be a vector of MONOTONIC entries\n');
                else
                    h_inv = 1./h;
                    % Construct a sparse tridiagonal matrix with the interval lengths.
                    % This matrix will be used to ontain the values of the second
                    % derivatives at the supplied points. Let's call it M
                    % We are trying to obtain the spline equation as:
                    % s''i(x) = (z_(i-1)/h_i)*(x_i-x) + (z_i/h_i)*(x-x_(i-1))  
                    % si(x) = (z_(i-1)/h_i)*(x_i-x)^3 + (z_i/h_i)*(x-x_(i-1))^3 + ...
                    %   Ci*(x_i-x) + Di*(x-x_(i-1))  
                    % Enforcing continuity of the second derivative at the boundary
                    % points gives us the valu
                    
                    M =  spdiags([h(2:end)' 2*(h(1:end-1)+h(2:end))' h(1:end-1)'], ...
                        -1:1, n_plus_1-2, n_plus_1-2);

                    % This is a very interesting way to create a triadiagonal matrix.
                    % M = spdiags(B, d, m, n) creates a sparse mxn matrix by taking the
                    % columns of B and placing them as diagonals (indexed by d). The
                    % -1:1 here means that the first column of B is placed on the
                    % digonal below the center diagonal, second column along the main
                    % diagonal and the third column along the diagonal above it.
                    
                    v_int = Obj.f_vals(Obj.colons{:}, 2:n_plus_1-1);    % MATLAB doesn't take in 'end' with cell-array
                    v_int_plus_1 = Obj.f_vals(Obj.colons{:}, 3:n_plus_1);
                    v_int_minus_1 = Obj.f_vals(Obj.colons{:}, 1:n_plus_1-2);  
                    % Shifted versions of the function values at the internal points; 
                    % We will solve Mz = b
                    
                    diff_dims = ndims(h) - ndims(Obj.f_vals);
                    h = shiftdim(h, diff_dims);
                    h_inv = shiftdim(h_inv, diff_dims);
                    b = 6.0 * ( times(v_int_minus_1, h_inv(Obj.colons{:}, 1:n_plus_1-2)) + ...
                        times(v_int_plus_1, h_inv(Obj.colons{:}, 2:n_plus_1-1)) - ...
                        times(v_int, (h_inv(Obj.colons{:}, 1:n_plus_1-2) + h_inv(Obj.colons{:}, 2:n_plus_1-1))) );
                    
                    % TODO: Can write my own code for LU decomposition of the tridiagonal
                    % matrix. It has a pretty nice structure (and is also diagonally
                    % dominant, so no pivots are required for solving crout-LU (CHECK)
            
                    % We would like to ideally use a sparse solver here because of the nice structure, but MATLAB keeps
                    % runninf into trouble with multi-dimensional matrices. We will use the inverse operation instead
                    % z = mrdivide(b, M);    % Generate a matrix of coeffcients

                    M_inv = full( inv(M) );
                    z = zeros( size(b) );

                    for c_index = 1 : size(M, 2)
                        z(Obj.colons{:}, c_index) = sum( times(b, shiftdim(M_inv(:, c_index)', diff_dims)), ndims(b) );
                    end

                    s_z = size(z);
                    nd_z = length(s_z);
                    z = cat(nd_z, zeros([s_z(1:end-1) 1]), z, zeros([s_z(1:end-1) 1])); 
                        % Natural spline assumes that the second derivative at the two
                        % end points is 0, which has been put in here

                    % Time to compute the coefficients nxCi and nxDi (as vectors C, D);
                    C__ = times( Obj.f_vals(Obj.colons{:}, 1:n_plus_1-1), h_inv ) - ...
                        times( z(Obj.colons{:}, 1:n_plus_1-1), h )/6.0;
                    D__ = times( Obj.f_vals(Obj.colons{:}, 2:n_plus_1  ), h_inv ) - ...
                        times( z(Obj.colons{:}, 2:n_plus_1), h )/6.0;

                    % RESHAPE z, C and D to be of the same size as the input matrix V
                    % (just to maintain consistency -- spline can be reapplied on the
                    % spline coefficients in this way)

                    zl = (1/6.0) * times( z(Obj.colons{:}, 1:n_plus_1-1), h_inv ); 
                    zr = (1/6.0) * times( z(Obj.colons{:}, 2:n_plus_1), h_inv );    

                    Obj.coeffs = zeros([Obj.op_dims 4 n_plus_1+1]);

                    % The formulation with zl, zr, C and D requires, as an input, (x -
                    % x_i) and (x_i+1 - x). This formulation does not extend so
                    % naturally for extrapolation, so we will move to local coordinates
                    % w.r.t the left element (left of the first element needs to be
                    % handled as a special case). The requires the coefficients to be
                    % multiplied by d^3, d^2, d, 1. The coefficients can be easily
                    % calculated if zl, zr, C and D are given.
                    Obj.coeffs(Obj.colons{:},1,2:n_plus_1) = zr - zl; % d^3
                    Obj.coeffs(Obj.colons{:},2,2:n_plus_1) = 3.0 * times(h, zl); % d^2
                    Obj.coeffs(Obj.colons{:},3,2:n_plus_1) = (D__ - C__) - 3.0 * times(h.^2, zl); % d
                    Obj.coeffs(Obj.colons{:},4,2:n_plus_1) = times(zl, h.^3) + times(h, C__); % 1

                    %% EXTRAPOLATION OFFSET IN coeffs.e(4,{1,end},colons{:})
                    Obj.coeffs(Obj.colons{:},4,1) =   zl(Obj.colons{:}, 1) * (h(1)^3) + ...
                                                C__(Obj.colons{:}, 1) * h(1);    % CnXn
                    Obj.coeffs(Obj.colons{:},4,n_plus_1+1) = zr(Obj.colons{:}, n_plus_1-1) * (h(end)^3) + ...
                                                D__(Obj.colons{:}, n_plus_1-1) * h(end);% DnXn

                    % EXTRAPOLATION SLOPE (RHS -- OK, LHS -- (hack) NEGATIVE SLOPE)
                    Obj.coeffs(Obj.colons{:},3,1) =   D__(Obj.colons{:}, 1) - C__(Obj.colons{:}, 1) + ...
                                                - 3.0 * zl(Obj.colons{:}, 1) * h(1);
                    Obj.coeffs(Obj.colons{:},3,n_plus_1+1) = D__(Obj.colons{:}, n_plus_1-1) - C__(Obj.colons{:}, n_plus_1-1) + ...
                                                + 3.0 * zr(Obj.colons{:}, n_plus_1-1) * h(end);
                end
            end    
        end

        function [val, der] = computeWithDer(Obj, x_in, ~)
            [x_in, dx_out] = Obj.i_pts.rescaleShiftInput(x_in);
            % We need to find the bins i: xq lies in (x_(i-1) - x_i)  
            %[~, index] = histc(x_in, [-Inf; Obj.getPtAt(1); Obj.i_pts.getPts(1); Inf]); 
            %if ( isa(Obj.i_pts, 'UniformPoints') )
            %    % Works for uniformly distributed points
            %    index = 1 + floor((x_in - Obj.i_pts.bounds(2,1))/Obj.i_pts.step_size(1));
            %    index = min( max(1, index), Obj.i_pts.n_pts(1)+1 );
            %else
                index = findInSorted(x_in, [Inf; Obj.i_pts.pts{1}; -Inf]);
            %end

            % Move to local coordinates, i.e., calculate the distances (and their cubes)
            % from the two end-points. Keep in mind that "index" is calculated by appending
            % Inf and -Inf to the set of sample points. Substract 1 to get the index of the
            % actual local coordinate

            dl_1 = x_in - Obj.i_pts.pts{1}(max(index-1,1));
            dl_2 = dl_1.*dl_1;
            dl_3 = dl_1.*dl_2;
            
            % Calculate the cubic function and the derivatives
            val =   dl_3*Obj.coeffs(Obj.colons{:},1,index) + ...
                    dl_2*Obj.coeffs(Obj.colons{:},2,index) + ...
                    dl_1*Obj.coeffs(Obj.colons{:},3,index) + ...
                         Obj.coeffs(Obj.colons{:},4,index);
            if (nargout > 1)
                der =   dx_out * ( ...
                        3*dl_2*Obj.coeffs(Obj.colons{:},1,index) + ...
                        2*dl_1*Obj.coeffs(Obj.colons{:},2,index) + ...
                               Obj.coeffs(Obj.colons{:},3,index));
            end
        end

        function der = secondDerivative(Obj, x_in)
            der = 0; % TODO
        end

        function der = firstDerivativeAtPt(Obj, pt_index)
            der = 0; % TODO
        end

        function der = secondDerivativeAtPt(Obj, pt_index)
            der = 0; % TODO
        end
    end
end