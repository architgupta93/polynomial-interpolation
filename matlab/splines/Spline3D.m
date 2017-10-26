classdef Spline3D < SplineInterpolant
% SPLINE3D
% Trial implementation for a 3D spline interpolant
% Interpolant will take in a vector of size 3 and evaluate a function of these
% three arguments.

    properties (Access = protected)
        
    end

    methods (Access = public)
        function Obj = Spline3D(varargin)
        % function Obj = Spline3D(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % Class Constructor
            Obj = Obj@SplineInterpolant(varargin{:});
            if ( isempty(varargin) )
                % Empty class instances are used by the save/load functions to
                % load the relevant class members into
                return;
            end

            % The parent classes have already done the bulk of the job in the call
            % Obj = Obj@SplineInterpolant(varargin{:}) for us.
            % Data has been sampled and stored in the field Obj.f_vals
            % The rightmost index in f_vals corresponds to the data for the
            % first input variable and as we move to the right, subsequence
            % variables are indexed.

            coeffs2D = zeros([Obj.op_dims, 4, 4, Obj.i_pts.getNPts(1)+1, ...
                Obj.i_pts.getNPts(2)+1, Obj.i_pts.getNPts(3)]);
            % Yes, that is a lot of coefficients if you are thinking along the
            % lines but what can one do

            % The following code is quite ugly but I would like this to work
            % before I write some beautifule code generalized for N-Dimensions

            % Let's deal with initial dimension by using Spline2D
            for di_index = 1:Obj.i_pts.getNPts(3)
                Si = Spline2D(Obj.f_vals(Obj.colons{:}, :, :, di_index), 2, ...
                    {Obj.bounds{1, 1}, Obj.bounds{1, 2}}, ...
                    Obj.order(1, 1:2), varargin{end});
                %       i_type_or_x_vals ^^

                % local spline object
                coeffs2D(Obj.colons{:},:,:,:,:,di_index) = Si.coeffs;
            end

            % We need to permute the coefficients a little since for the
            % righmost dimension, we will get size [ ..., 4, n_pts], whereas,
            % at the end of it all, we would like to put all the 4s together
            % (if that makes any sense to you reader)

            % Keep the output intact, and swap
            % [4_n2, 4_n1, n1, n2, 4_n3, n3] to [4_n3, 4_n2, 4_n1, n1, n2, n3]...
            % The coefficients go: [5 1 2 3 4 6]
            perm_order = [ [1 : g_dims(Obj.op_dims)], [g_dims(Obj.op_dims) + [5 1 2 3 4 6]] ];
            Sij = Spline1D(coeffs2D, 1, {Obj.bounds{1,3}}, Obj.order(1,3), varargin{end});
            Obj.coeffs = permute(Sij.coeffs, perm_order);
        end

        function [val, der] = computeWithDer(Obj, x_in, ~)
            val = 0;
            der = 0;
        end
        
    % end methods
    end
    
% end classdef
end
