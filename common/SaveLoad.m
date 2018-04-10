classdef SaveLoad < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classdef SaveLoad < handle
% Class for manually implementing Save/Load interface for large objects.
% MATLAB provides a very flexible interface for saving/loading objects,
% especially instances of user-defined classes which is quite non-trivial.
%
%           However, Octave doesn't extend this functionality.
%
% Our implementation of Save/Load lets us extend the functionality of
% saving/loading class instance by simply inheriting this class.
%
% REQUIREMENT: If one or more of the members of a child class do not natively
% support Save/Load functions then they should also be derived from this
%
% Author: Archit Gupta
% Date: October 03, 2017
%
% SaveLoad class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Constant, Access =  protected)
        % Default matlab options
        load_opts = {'-mat'};

        % Having -v7 opyion in MATLAB lets us store objects that are larger
        % than 2 GB, or 8 maybe, I am not sure.
        save_opts = {'-v7'};
    end
end
