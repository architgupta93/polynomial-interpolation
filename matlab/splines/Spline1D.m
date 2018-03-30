classdef Spline1D < SplineInterpolant
    methods (Access = public)
        function Obj = Spline1D(varargin)
        % function Obj = Spline1D(f_vals, n_in_dims, bounds, order, i_type_or_x_vals)        % Class Constructor
            Obj = Obj@SplineInterpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            Obj.setupCoefficients();
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
            % In MATLAB, all singleton dimensions to the right are ignored, so
            % in the expression below, the dimension corresponding to 'index'
            % is ignored.
            ext_dims = length(Obj.colons);
            pj_vals  = shiftdim([dl_3, dl_2, dl_1, 1], 1-ext_dims);
            val      = sum(pj_vals .* ...
                       Obj.coeffs(Obj.colons{:},:,index), 1+ext_dims);

            if (nargout > 1)
                dpj_vals = shiftdim([3*dl_2, 2*dl_1, 1], 1-ext_dims);
                der      = dx_out * sum(dpj_vals .* ...
                           Obj.coeffs(Obj.colons{:},1:3,index), 1+ext_dims);
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
