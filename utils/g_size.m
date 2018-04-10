function n__sz = g_size(in_x)
    if (isvector(in_x))
        n__sz = length(in_x);    
    else
        n__sz = size(in_x);
    end
end
