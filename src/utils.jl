pad_to_ndim(vec) = begin
    out = @MVector zeros(vofi_real, NDIM)
    for i in 1:min(length(vec), NDIM)
        out[i] = vofi_real(vec[i])
    end
    out
end
