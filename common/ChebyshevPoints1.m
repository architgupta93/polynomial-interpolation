classdef ChebyshevPoints1 < ChebyshevPoints
    methods (Access = protected)

    end

    methods (Access = public)
        % Class constructor
        function Obj = ChebyshevPoints1(n_in_dims, order, bounds, roots_or_extremas)
            Obj = Obj@ChebyshevPoints(n_in_dims, order, bounds);
            if (nargin < 4) % If nothing has been supplied, use extremas for
                            % Chebyshev points of the first kind
                Obj.roots_or_extremas = false;
            else
                Obj.roots_or_extremas = roots_or_extremas;
            end

            pts = Obj.generate(Obj.order);
            Obj.recenterAndSave(pts);
            Obj.cheb_type = 1;
        end

        function [pts, n_pts] = generate(Obj, order)
            pts = cell(Obj.n_in_dims, 1);
            if (Obj.roots_or_extremas)  % Roots of T_(n_pts+1)
                                        % Cannot be used recursively
                % Solving for cos( (n_pts+1) acos(x) ) = 0
                n_pts = 2.^order;
                % Initialize the chebyshev points of type 1
                for d_i = 1:Obj.n_in_dims
                    pspace = (0.5/n_pts) + linspace( ...
                        0, (n_pts(d_i)-1)/n_pts(d_i), n_pts(d_i))';
                    pts{d_i} = cos( pi * pspace );
                end
            else    % Extremas of T_(n_pts)
                % Solving for | cos( (n_pts) * acos(x) ) | = 1
                % acos(x) = k * pi / (n_pts)
                n_pts = 2 .^ order + 1;
                for d_i = 1:Obj.n_in_dims
                    if (n_pts(d_i) <= 1)
                        pspace = 0.5;
                    else
                        pspace = linspace( 0, 1, n_pts(d_i) )';
                    end
                    pts{d_i} = cos( pi * pspace );
                end
            end
        end

        function wts = getBLIWeights(Obj, scale)
            % See UniformPoints.m for a definition of BLI points. For the case
            % of Chebyshev points of type 2, the expression for wts is a little
            % too trivial (at least in one dimension):
            %
            %           wj = (-1)^j * 1/2 if j is an endpoint
            %           wj = (-1)^j * 1 otherwise
            %
            % TODO: Note that the implementation of getBLIWeights for both
            % UniformPoints and ChebyshevPoints might be incorrect for 2 or
            % higher dimensions (Actually, I am certain that it is wrong). Will
            % look into it when the time comes for implementing BLI in higher
            % dimensions

            wts = cell(Obj.n_in_dims, 1);
            for d_i = 1:Obj.n_in_dims
                wts{d_i} = 1 - 2*mod( [0:Obj.n_pts(d_i)-1]', 2);
                wts{d_i}(1) = 0.5;
                wts{d_i}(end) = wts{d_i}(end)/2;
                if (nargin > 1)
                    wts{d_i} = wts{d_i} * scale(d_i);
                end
            end
        end
    end
end
