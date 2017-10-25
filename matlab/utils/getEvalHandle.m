function f_handle = getEvalHandle(Obj)
%function f_handle = getEvalHandle(Obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ( isa(Obj, 'function_handle') )
        f_handle = Obj;
    elseif ( isa(Obj, 'Interpolant') | isa(Obj, 'PiecewiseInterpolant') )
        f_handle = @(x_in) Obj.computeWithDer( x_in );
    else
        error('ERROR: Invalid object type for extracting function handle');
    end
end
