function M = flipdims(M)

    for k = 1:ndims(M)
        M = flipdim(M,k);
    end
end