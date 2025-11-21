function vofi_check_boundary_line(impl_func, par, x0, h0, f0, xfs_pt, n0)
    nx = ones(vofi_int, NSE)
    ny = ones(vofi_int, NSE)
    sidedirx = (1.0, 0.0, 0.0)
    sidediry = (0.0, 1.0, 0.0)
    fse = zeros(vofi_real, NSE)
    x1 = zeros(vofi_real, NDIM)
    xfsl = MinData()
    check_dir = -1

    for i in 0:1
        for j in 0:1
            if n0[i + 1, j + 1] > 0
                if ny[i + 1] > 0
                    ny[i + 1] = 0
                    fse[1] = f0[i + 1, 1]
                    fse[2] = f0[i + 1, 2]
                    for k in 1:NDIM
                        x1[k] = x0[k] + i * sidedirx[k] * h0[1]
                    end
                    consi = vofi_check_side_consistency(impl_func, par, x1, sidediry, fse, h0[2])
                    if consi != 0
                        f2pos = consi
                        sign_change = vofi_get_segment_min(impl_func, par, x1, sidediry, fse, xfsl, h0[2], f2pos)
                        if sign_change != 0
                            copy!(xfs_pt, xfsl)
                            xfs_pt.isc[1] = 1
                            xfs_pt.isc[i + 2] = 1
                            check_dir = 0
                        end
                    end
                end
                if nx[j + 1] > 0
                    nx[j + 1] = 0
                    fse[1] = f0[1, j + 1]
                    fse[2] = f0[2, j + 1]
                    for k in 1:NDIM
                        x1[k] = x0[k] + j * sidediry[k] * h0[2]
                    end
                    consi = vofi_check_side_consistency(impl_func, par, x1, sidedirx, fse, h0[1])
                    if consi != 0
                        f2pos = consi
                        sign_change = vofi_get_segment_min(impl_func, par, x1, sidedirx, fse, xfsl, h0[1], f2pos)
                        if sign_change != 0
                            copy!(xfs_pt, xfsl)
                            xfs_pt.isc[1] = 1
                            xfs_pt.isc[j + 2] = 1
                            check_dir = 1
                        end
                    end
                end
                n0[i + 1, j + 1] = 0
            end
        end
    end

    return check_dir
end

function vofi_check_boundary_surface(impl_func, par, x0, h0, f0, xfs, n0)
    nx = ones(vofi_int, NSE)
    ny = ones(vofi_int, NSE)
    nz = ones(vofi_int, NSE)
    sidedirx = (1.0, 0.0, 0.0)
    sidediry = (0.0, 1.0, 0.0)
    sidedirz = (0.0, 0.0, 1.0)
    fve = zeros(vofi_real, NVER)
    x1 = zeros(vofi_real, NDIM)
    xfsl = MinData()
    check_dir = -1

    for i in 0:1
        for j in 0:1
            for k in 0:1
                if n0[i + 1, j + 1, k + 1] > 0
                    if nx[i + 1] > 0
                        nx[i + 1] = 0
                        fve[1] = f0[i + 1, 1, 1]
                        fve[2] = f0[i + 1, 2, 1]
                        fve[3] = f0[i + 1, 1, 2]
                        fve[4] = f0[i + 1, 2, 2]
                        for m in 1:NDIM
                            x1[m] = x0[m] + i * sidedirx[m] * h0[1]
                        end
                        ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sidediry, sidedirz, fve)
                        if ipsc.consi != 0
                            sign_change = vofi_get_face_min(impl_func, par, x1, h0, sidediry, sidedirz, fve, xfsl, ipsc)
                            if sign_change != 0
                                copy!(xfs[1], xfsl)
                                xfs[1].isc[1] = 1
                                xfs[1].isc[i + 2] = 1
                                check_dir = 0
                            end
                        end
                    end
                    if ny[j + 1] > 0
                        ny[j + 1] = 0
                        fve[1] = f0[1, j + 1, 1]
                        fve[2] = f0[2, j + 1, 1]
                        fve[3] = f0[1, j + 1, 2]
                        fve[4] = f0[2, j + 1, 2]
                        for m in 1:NDIM
                            x1[m] = x0[m] + j * sidediry[m] * h0[2]
                        end
                        ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sidedirx, sidedirz, fve)
                        if ipsc.consi != 0
                            sign_change = vofi_get_face_min(impl_func, par, x1, h0, sidedirx, sidedirz, fve, xfsl, ipsc)
                            if sign_change != 0
                                copy!(xfs[2], xfsl)
                                xfs[2].isc[1] = 1
                                xfs[2].isc[j + 2] = 1
                                check_dir = 0
                            end
                        end
                    end
                    if nz[k + 1] > 0
                        nz[k + 1] = 0
                        fve[1] = f0[1, 1, k + 1]
                        fve[2] = f0[2, 1, k + 1]
                        fve[3] = f0[1, 2, k + 1]
                        fve[4] = f0[2, 2, k + 1]
                        for m in 1:NDIM
                            x1[m] = x0[m] + k * sidedirz[m] * h0[3]
                        end
                        ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sidedirx, sidediry, fve)
                        if ipsc.consi != 0
                            sign_change = vofi_get_face_min(impl_func, par, x1, h0, sidedirx, sidediry, fve, xfsl, ipsc)
                            if sign_change != 0
                                copy!(xfs[3], xfsl)
                                xfs[3].isc[1] = 1
                                xfs[3].isc[k + 2] = 1
                                check_dir = 0
                            end
                        end
                    end
                    n0[i + 1, j + 1, k + 1] = 0
                end
            end
        end
    end

    return check_dir
end

function vofi_check_secondary_side(impl_func, par, x0, h0, pdir, sdir, f0, xfs_pt, fth)
    hs = zero(vofi_real)
    for i in 1:NDIM
        hs += sdir[i] * h0[i]
    end
    x1 = zeros(vofi_real, NDIM)
    fse = zeros(vofi_real, NSE)
    xfsl = MinData()

    for k in 0:1
        fse[1] = f0[k + 1, 1]
        fse[2] = f0[k + 1, 2]
        if fse[1] * fse[2] < 0
            xfs_pt.isc[1] = 1
            xfs_pt.isc[k + 2] = -1
        else
            if abs(fse[1]) > fth && abs(fse[2]) > fth
                continue
            end
            for i in 1:NDIM
                x1[i] = x0[i] + k * pdir[i] * h0[i]
            end
            consi = vofi_check_side_consistency(impl_func, par, x1, sdir, fse, hs)
            if consi != 0
                f2pos = consi
                sign_change = vofi_get_segment_min(impl_func, par, x1, sdir, fse, xfsl, hs, f2pos)
                if sign_change != 0
                    copy!(xfs_pt, xfsl)
                    xfs_pt.isc[1] = 1
                    xfs_pt.isc[k + 2] = 1
                end
            end
        end
    end
end

function vofi_check_secter_face(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfs_pt, fth)
    xfs_pt.isc .= 0
    x1 = zeros(vofi_real, NDIM)
    fve = zeros(vofi_real, NVER)
    xfsl = MinData()

    for m in 0:1
        np0 = nm0 = 0
        fve[1] = f0[m + 1, 1, 1]; np0 += fve[1] > 0; nm0 += fve[1] < 0
        fve[2] = f0[m + 1, 2, 1]; np0 += fve[2] > 0; nm0 += fve[2] < 0
        fve[3] = f0[m + 1, 1, 2]; np0 += fve[3] > 0; nm0 += fve[3] < 0
        fve[4] = f0[m + 1, 2, 2]; np0 += fve[4] > 0; nm0 += fve[4] < 0

        if nm0 * np0 == 0
            if !(abs(fve[1]) > fth && abs(fve[2]) > fth && abs(fve[3]) > fth && abs(fve[4]) > fth)
                for i in 1:NDIM
                    x1[i] = x0[i] + m * pdir[i] * h0[i]
                end
                ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sdir, tdir, fve)
                if ipsc.consi != 0
                    sign_change = vofi_get_face_min(impl_func, par, x1, h0, sdir, tdir, fve, xfsl, ipsc)
                    if sign_change != 0
                        copy!(xfs_pt, xfsl)
                        xfs_pt.isc[1] = 1
                        xfs_pt.isc[m + 1] = 1
                    end
                end
            end
        end
    end
    return nothing
end

function vofi_check_tertiary_side(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfs, fth)
    ht = sum(tdir .* h0)
    for idx in 1:4
        xfs[idx].isc .= 0
    end
    x1 = zeros(vofi_real, NDIM)
    fse = zeros(vofi_real, NSE)

    for m in 0:1
        for n in 0:1
            fse[1] = f0[m + 1, n + 1, 1]
            fse[2] = f0[m + 1, n + 1, 2]
            l0 = 2 * m + n + 1
            if fse[1] * fse[2] < 0
                xfs[l0].isc[1] = 1
                xfs[l0].isc[2] = -1
            else
                if !(abs(fse[1]) > fth && abs(fse[2]) > fth)
                    for i in 1:NDIM
                        x1[i] = x0[i] + m * pdir[i] * h0[i] + n * sdir[i] * h0[i]
                    end
                    consi = vofi_check_side_consistency(impl_func, par, x1, tdir, fse, ht)
                    if consi != 0
                        xfsl = MinData()
                        f2pos = consi
                        sign_change = vofi_get_segment_min(impl_func, par, x1, tdir, fse, xfsl, ht, f2pos)
                        if sign_change != 0
                            xfs[l0].xval .= xfsl.xval
                            xfs[l0].sval = xfsl.sval
                            xfs[l0].fval = xfsl.fval
                            xfs[l0].isc[1] = 1
                            xfs[l0].isc[2] = 1
                        end
                    end
                end
            end
        end
    end
    return nothing
end
