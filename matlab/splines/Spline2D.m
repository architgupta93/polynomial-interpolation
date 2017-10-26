 classdef Spline2D < SplineInterpolant
    methods (Access = public)
        function Obj = Spline2D(varargin)
        % function Obj = Spline2D(f_vals, in_dims, bounds, order, i_type_or_x_vals)
        % Class Constructor
            Obj = Obj@SplineInterpolant(varargin{:});
            if ( isempty(varargin) )
                return;
            end

            % Since the first entry in vij(;,:,:,...) is indexed by xi, we start by
            % constructing splines with xi as a variable (different fixed xjs)

            coeffs1D = zeros([Obj.op_dims 4 Obj.i_pts.getNPts(1)+1 Obj.i_pts.getNPts(2)]);
            for jIndex = 1:Obj.i_pts.getNPts(2)
                Sj = Spline1D(Obj.f_vals(Obj.colons{:}, :, jIndex), 1, Obj.bounds{1, 1}, Obj.order(1, 1), varargin{end});
                %                                                                           i_type_or_x_vals ^^
                % local spline object
                coeffs1D(Obj.colons{:},:,:,jIndex) = Sj.coeffs;
            end

            % Now we will flip the coefficient tensor around so that we can do a
            % spline interpolation on the coefficients

            % Swap the last two entries in this vector
            % perm_order = [1 : ndims(coeffs1D)];
            % perm_order(end) = perm_order(end-1); 
            % perm_order(end-1) = ndims(coeffs1D);
            % coeffs1D = permute(coeffs1D, perm_order);

            % Spline1D returns size [Obj.op_dims 4 1+Obj.i_pts.getNPts(1) 4 1+Obj.i_pts.getNPts(2)].
            % We need to  swap this to [Obj.op_dims 4 4 1+Obj.i_pts.getNPts(1) 1+Obj.i_pts.getNPts(2)]
            % Where the 2 4s have also been swapped, i.e., from [1 2 3 4], we go to [3 1 2 4]
            perm_order = [ [1 : g_dims(Obj.op_dims)], [g_dims(Obj.op_dims) + [3 1 2 4]] ];
            Obj.coeffs = zeros([Obj.op_dims 4 4 1+Obj.i_pts.getNPts(1) 1+Obj.i_pts.getNPts(2)]);
            for iIndex = 1:Obj.i_pts.n_pts(1)+1
                coeffs_xi_iIndex = coeffs1D(Obj.colons{:},:,iIndex,:);

                % All the coefficients corresponding to xi
                Sij = Spline1D(coeffs_xi_iIndex, 1,  Obj.bounds{1, 2}, Obj.order(1, 2), varargin{end});
                %                                                           i_type_or_x_vals ^^
                Obj.coeffs(Obj.colons{:},:,:,iIndex,:) = permute(Sij.coeffs, perm_order);
            end

            % Notice that spline_1d implementation already provides the correct
            % coefficients for slopes even for the case of extrapolation. The only
            % thing that is missing is the extrapolation bias values, which we can
            % supply easily

            % Lets turn the plane around and look at extrapolation beyond the
            % sampled xj boundary, controlled by coeffSpline.v([zl,zr,C,D], 1, ...
            % :,:,colons{:}). We apply the same strategy as before and let the
            % bicubic coefficients grow linearly in xj. However, since  these
            % coeffecients have been obtained by spline_1d, they have already been
            % organized correctly
        end

        function [val, der] = computeWithDer(Obj, x_in, ~)
            [x_in, dx_out] = Obj.i_pts.rescaleShiftInput(x_in);
            %if ( isa(Obj.i_pts, 'UniformPoints') )
            %    % Works for uniformly distributed points
            %    i = 1 + floor((x_in(1) - Obj.i_pts.bounds(2,1))/Obj.i_pts.step_size(1));
            %    j = 1 + floor((x_in(2) - Obj.i_pts.bounds(2,2))/Obj.i_pts.step_size(2));
            %    i = min( max(1, i), Obj.i_pts.n_pts(1)+1 );
            %    j = min( max(1, j), Obj.i_pts.n_pts(2)+1 );
            %else
                i = findInSorted(x_in(1), [Inf; Obj.i_pts.pts{1}; -Inf]);
                j = findInSorted(x_in(2), [Inf; Obj.i_pts.pts{2}; -Inf]);
            %end
            
            dil_1 = x_in(1) - Obj.i_pts.pts{1}(max(i-1,1));
            djl_1 = x_in(2) - Obj.i_pts.pts{2}(max(j-1,1));
                                            % Distance of our point in the local
                                            % coordinates, i.v., from the left end-point
                                            % of the segments that it lies in
            dil_2 = dil_1.*dil_1;
            dil_3 = dil_1.*dil_2;

            djl_2 = djl_1.*djl_1;
            djl_3 = djl_1.*djl_2;
                    

            % The coefficients are a OutputSizex4x(xi)x4x(xj) tensor. The 4x4 matrix
            % essentially provides the 16 coefficients required for computing the
            % function value in the rectangle (or extrapolation region) given by
            % (xi(i), xi(j)) and (xi(i+1),xj(j+1)) 

            cubicCoeffs3 =  djl_3*Obj.coeffs(Obj.colons{:},1,1,i,j)+ ...
                            djl_2*Obj.coeffs(Obj.colons{:},2,1,i,j)+ ...        
                            djl_1*Obj.coeffs(Obj.colons{:},3,1,i,j)+ ...        
                                  Obj.coeffs(Obj.colons{:},4,1,i,j);
                                                                  
            cubicCoeffs2 =  djl_3*Obj.coeffs(Obj.colons{:},1,2,i,j)+ ...
                            djl_2*Obj.coeffs(Obj.colons{:},2,2,i,j)+ ...        
                            djl_1*Obj.coeffs(Obj.colons{:},3,2,i,j)+ ...        
                                  Obj.coeffs(Obj.colons{:},4,2,i,j);
                                                                  
            cubicCoeffs1 =  djl_3*Obj.coeffs(Obj.colons{:},1,3,i,j)+ ...
                            djl_2*Obj.coeffs(Obj.colons{:},2,3,i,j)+ ...        
                            djl_1*Obj.coeffs(Obj.colons{:},3,3,i,j)+ ...        
                                  Obj.coeffs(Obj.colons{:},4,3,i,j);
                                                                  
            cubicCoeffs0 =  djl_3*Obj.coeffs(Obj.colons{:},1,4,i,j)+ ...
                            djl_2*Obj.coeffs(Obj.colons{:},2,4,i,j)+ ...        
                            djl_1*Obj.coeffs(Obj.colons{:},3,4,i,j)+ ...        
                                  Obj.coeffs(Obj.colons{:},4,4,i,j);

            val =   cubicCoeffs3*dil_3 + ...
                    cubicCoeffs2*dil_2 + ...
                    cubicCoeffs1*dil_1 + ...
                    cubicCoeffs0;

            if (nargout > 1)
                der = zeros([Obj.op_dims Obj.in_dims]);
                der(Obj.colons{:}, 1) = dx_out(1) * ( ...
                                        3*cubicCoeffs3*dil_2 + ...
                                        2*cubicCoeffs2*dil_1 + ...
                                        cubicCoeffs1);

                dcubicCoeffs3 = ...
                    3*djl_2*Obj.coeffs(Obj.colons{:},1,1,i,j)+ ...
                    2*djl_1*Obj.coeffs(Obj.colons{:},2,1,i,j)+ ...
                            Obj.coeffs(Obj.colons{:},3,1,i,j);
                                                            
                dcubicCoeffs2 = ...                         
                    3*djl_2*Obj.coeffs(Obj.colons{:},1,2,i,j)+ ...
                    2*djl_1*Obj.coeffs(Obj.colons{:},2,2,i,j)+ ...
                            Obj.coeffs(Obj.colons{:},3,2,i,j);
                                                            
                dcubicCoeffs1 = ...                         
                    3*djl_2*Obj.coeffs(Obj.colons{:},1,3,i,j)+ ...
                    2*djl_1*Obj.coeffs(Obj.colons{:},2,3,i,j)+ ...
                            Obj.coeffs(Obj.colons{:},3,3,i,j);
                                                            
                dcubicCoeffs0 = ...                         
                    3*djl_2*Obj.coeffs(Obj.colons{:},1,4,i,j)+ ...
                    2*djl_1*Obj.coeffs(Obj.colons{:},2,4,i,j)+ ...
                            Obj.coeffs(Obj.colons{:},3,4,i,j);

                der(Obj.colons{:}, 2) = dx_out(2) * ( ... 
                                 dcubicCoeffs3*dil_3 + ...
                                 dcubicCoeffs2*dil_2 + ...
                                 dcubicCoeffs1*dil_1 + ...
                                 dcubicCoeffs0);
            end
        end

        function der = secondDerivative(Obj, x_in)
            der = 0; % TODO
        end

        function der = firstDerivativeAtPt(Obj, x_in)
            der = 0; % TODO
        end

        function der = secondDerivativeAtPt(Obj, x_in)
            der = 0; % TODO
        end
    end
end
