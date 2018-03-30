classdef ChebSeries < Interpolant
% The ChebSeries class is used to represent a Chebyshev series. It represents
% the following polynomial:
%       y = a0 * C0(x) + a1 * C1(x) + ... + an * Cn(x),
% Where, Ci(x) is the chebyshev polynomial of some particular type.
% For example, for Chebyshev polynomials of the first kind, the polynomial Ci(x)
% is given by Ci(x) = cos( i * acos(x) )
%
% Author: Archit Gupta (Feb 06, 2017)

    properties (Constant)
        ERR_TH = 1e-14; % Error threshold to decide the number of series
                        % coefficients
    end

    properties (Access = protected)
        cheb_type = 0;
    end

    methods (Access = public)
        function Obj = ChebSeries(varargin)
        % function Obj = ChebSeries(f_vals, in_dims, bounds, order, i_type_or_x_vals)

            % Call the parent class constructor 
            varargin{5} = 'chebyshev';
            Obj         = Obj@Interpolant(varargin{:});

            % TODO: Make sure that we have an instance of ChebyshevPoints
            % analyst class in i_pts otherwise getting SeriesCoefficients
            % doesn't make any sense

            % Get Chebyshev series coefficients
            Obj.coeffs = Obj.i_pts.getSeriesCoeffs(Obj.f_vals');
        end

        function pts = generate(Obj, order)
            pts = Obj.cheb_pts.generate(order);
        end

        function plot(Obj)
            figure();
            if (Obj.i_pts.in_dims == 1)
                semilogy(abs(Obj.coeffs), 'LineWidth', 2.0);
                xlabel('Coefficient Index (n)');
                ylabel('Coefficient Value (log scale)');
            elseif (Obj.i_pts.in_dims == 2)
                surf( abs(Obj.coeffs) + 1e-16 );
                xlabel('Coeff Index (x)');
                ylabel('Coeff Index (y)');
                zlabel('Coeff Value (log scale)');
                set(gca, 'ZScale', 'log');
                colormap Autumn;
            else
                fprintf(['Cannot plot higher that 2-dimensional data\n']);
            end
            set(gca, 'FontSize', 28);
            grid on;
        end

        function [val, der] = computeWithDer(Obj, x_in)
            y_out = Obj.i_pts.evaluatePolynomial(x_in);    
            val = bsxfun( @times, Obj.coeffs, tensorProduct(y_out) );
            for d_i = 1:Obj.in_dims
                val = squeeze( sum(val) );
            end

            % TODO: Compute the derivative too!
            der = 0;
        end
    end
end
