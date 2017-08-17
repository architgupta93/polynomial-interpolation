classdef SplineInterpolant < Interpolant
% SplineInterpolant: Base class for various spline interpolants, currently including Spline1D and Spline2D.
% This class provides class member called "coeffs" which is used for storing Spline interpolants in both the classes.
% The class also provides a common interface for defining Save and Load methods for the class
    properties (GetAccess = public, SetAccess = protected)
        coeffs = [];    
    end

    methods (Access = public)
        function Obj = SplineInterpolant(varargin)
            % HACK HACK HACK: We use PiecewiseBLI at the same level as Spline2D. Now spline, being an interpolant, takes
            % in a matrix as bounds (as it can only be one piece at the moment. However, since it is used in the same
            % API as PiecewiseBLI, a cell array is passed in at this level.
            if ( length(varargin) > 0 ) % Then everything HAS to be present
                if ( isa(varargin{3}, 'cell') )
                    varargin{3} = cell2mat(varargin{3});
                end
            end
            Obj = Obj@Interpolant(varargin{:});
        end

        function load(Obj, filename)
            % Call load on the superclass
            load@Interpolant(Obj, filename);

            % Load the remaining class
            % DEBUG
            % fprintf(2, 'Loading Spline Add-ons...\n');
            % END DEBUG
            load([filename, '.sco'], ... Spline COefficients
                'coeffs', ...
                Obj.load_opts{:});
            Obj.coeffs = coeffs;
        end

        function save(Obj, filename)
            % Call save on the superclass object
            save@Interpolant(Obj, filename);

            % Save the remaining class members as Add-on
            % DEBUG
            % fprintf(2, 'Saving Spline Add-ons...\n');
            % END DEBUG
            coeffs = Obj.coeffs;
            save([filename, '.sco'], ... Spline COefficients
                'coeffs', ...
                Obj.save_opts{:});
        end
    end
end
