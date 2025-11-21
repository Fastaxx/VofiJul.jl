pad_to_ndim(vec) = begin
    out = zeros(vofi_real, NDIM)
    vals = collect(vec)
    for i in 1:min(length(vals), NDIM)
        out[i] = vofi_real(vals[i])
    end
    out
end
