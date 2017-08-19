function n__dims = g_dims(in_x)
    n__dims = isvector(in_x) + ( ~isvector(in_x) * ndims(in_x) );
end
