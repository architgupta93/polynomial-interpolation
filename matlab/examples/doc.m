%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% (Piecewise) Polynomial Interpolation                                         %
% MATLAB/Octave Implementation                                                 %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some of the popular polynomial interpolation methods have been put together
% with a single access interface to make it easier to test different kinds of
% interpolants on a single problem. Besides interpolating the actual function
% values, this can also provide the derivative of the interpolant that is
% approximating the function (or values) that the user has supplied. This can
% be useful in several cases where the value for a function is known at a
% discrete set of points and an approximation of the derivative is required
% (say for finding zeros of the function). Presently, this includes:
%
% 1. Cubic splines
% 2. Lagrange Interpolant
%   - Piecewise Lagrange Interpolant
% 3. Barycentric Lagrange Interpolant
%   - Piecewise BLI
% 4. Discrete Cosine basis
% 
% Examples can be accessed in the example directory
%
% For an illustrative example on the usage of the package, run
% TODO
% >> help illustrative_example
%
% This will guide you through a simplified way to use the interpolation package
% on a function handle which takes a single vector (scalar works too) as an
% input and provides a multi-dimensional matrix as the output (can be a scalar,
% or vector too)
%
% For a detailed example on the various options that the top level functions
% take for interpolating function handles (same constraints on the function
% handle as described in illustrative_example), run
% TODO
% >> help detailed_example
%
% Both detailed_example and illustrative_example use splines to interpolate the
% supplied data. You can look around the examples to see the documentation for
% other interpolants as well
% TODO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Author: Archit Gupta
% Date: August 16, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
