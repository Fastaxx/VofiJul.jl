function vofi_order_dirs_1D(impl_func, par, x0, h0, f0)
    np0 = 0
    nm0 = 0
    icc = -1
    nmax0 = 2
    x1 = @MVector zeros(vofi_real, NDIM)
    MIN_GRAD = 1.0e-4

    x1[2] = x0[2]
    x1[3] = x0[3]
    
    # Evaluate at the two endpoints
    for i in 0:1
        x1[1] = x0[1] + i * h0[1]
        val = call_integrand(impl_func, par, x1)
        f0[i + 1] = val
        if val > 0
            np0 += 1
        elseif val < 0
            nm0 += 1
        end
    end

    # Compute simple gradient
    fgrad = (f0[2] - f0[1]) / h0[1]
    fgradmod = max(abs(fgrad), MIN_GRAD)
    hm = 0.5 * h0[1]
    fth = fgradmod * hm

    # Check if all endpoints have same sign
    if np0 * nm0 == 0
        np0 = nm0 = 0
        for i in 0:1
            f0mod = abs(f0[i + 1])
            if f0mod > fth
                if f0[i + 1] < 0
                    nm0 += 1
                else
                    np0 += 1
                end
            end
        end
        if nm0 == nmax0
            return 1  # fully inside
        elseif np0 == nmax0
            return 0  # fully outside
        end
        # Near boundary
        return nm0 > 0 ? 1 : 0
    end

    return -1  # interface crosses the cell
end

function vofi_order_dirs_2D(impl_func, par, x0, h0, pdir, sdir, f0, xfs_pt)
    n0 = @MMatrix zeros(Int, NSE, NSE)
    fc = @MMatrix zeros(vofi_real, NDIM, NDIM)
    hh = 0.5 .* h0
    np0 = 0
    nm0 = 0
    icc = -1
    check_dir = -1
    nmax0 = 4
    x1 = @MVector zeros(vofi_real, NDIM)
    fgrad = @MVector zeros(vofi_real, NSE)
    MIN_GRAD = 1.0e-4

    x1[3] = x0[3]
    for i in 0:1
        for j in 0:1
            x1[1] = x0[1] + i * h0[1]
            x1[2] = x0[2] + j * h0[2]
            val = call_integrand(impl_func, par, x1)
            f0[i + 1, j + 1] = val
            if val > 0
                np0 += 1
            elseif val < 0
                nm0 += 1
            end
        end
    end

    fgrad[1] = 0.5 * ((f0[2, 2] + f0[2, 1]) - (f0[1, 2] + f0[1, 1])) / h0[1]
    fgrad[2] = 0.5 * ((f0[2, 2] + f0[1, 2]) - (f0[2, 1] + f0[1, 1])) / h0[2]
    fgradsq = Sq2(fgrad)
    fgradmod = max(sqrt(fgradsq), MIN_GRAD)
    hm = max(hh[1], hh[2])
    fth = fgradmod * hm

    if np0 * nm0 == 0
        np0 = nm0 = 0
        for i in 0:1
            for j in 0:1
                f0mod = abs(f0[i + 1, j + 1])
                if f0mod > fth
                    n0[i + 1, j + 1] = 0
                    if f0[i + 1, j + 1] < 0
                        nm0 += 1
                    else
                        np0 += 1
                    end
                else
                    n0[i + 1, j + 1] = 1
                end
            end
        end
        if nm0 == nmax0
            return 1
        elseif np0 == nmax0
            return 0
        end
        check_dir = vofi_check_boundary_line(impl_func, par, x0, h0, f0, xfs_pt, n0)
        if check_dir < 0
            return nm0 > 0 ? 1 : 0
        end
    end

    have = 0.5 * (h0[1] + h0[2])
    for i in 0:1
        for j in 0:1
            fc[2 * i + 1, 2 * j + 1] = f0[i + 1, j + 1]
        end
    end
    for i in 1:NDIM
        x1[i] = x0[i] + hh[i]
    end
    fc[2, 2] = call_integrand(impl_func, par, x1)
    for i in 0:2:2
        x1[1] = x0[1] + i * hh[1]
        fc[i + 1, 2] = call_integrand(impl_func, par, x1)
    end
    x1[1] = x0[1] + hh[1]
    for j in 0:2:2
        x1[2] = x0[2] + j * hh[2]
        fc[2, j + 1] = call_integrand(impl_func, par, x1)
    end

    nc = @MMatrix zeros(Int, NDIM, NDIM)
    for i in 1:NDIM
        for j in 1:NDIM
            val = fc[i, j]
            nc[i, j] = val > 0 ? 1 : val < 0 ? -1 : 0
        end
    end

    nix = niy = 0
    for i in 0:2:2
        if nc[i + 1, 3] * nc[i + 1, 1] < 0
            niy += 1
        end
    end
    for j in 0:2:2
        if nc[3, j + 1] * nc[1, j + 1] < 0
            nix += 1
        end
    end

    fill!(pdir, 0)
    fill!(sdir, 0)
    if niy > nix
        jp, js = 2, 1
    elseif nix > niy
        jp, js = 1, 2
    else
        fgrad .= 0
        rwgt = (100.0, 50.0, 10.0, 2.0, 1.0)
        for i in 1:2
            for j in 1:2
                fx = 0.5 * ((fc[i + 1, j + 1] + fc[i + 1, j]) -
                            (fc[i, j + 1] + fc[i, j])) / hh[1]
                fy = 0.5 * ((fc[i + 1, j + 1] + fc[i, j + 1]) -
                            (fc[i + 1, j] + fc[i, j])) / hh[2]
                tmp = sqrt(fx^2 + fy^2)
                fx /= tmp
                fy /= tmp
                iwgt = abs(nc[i + 1, j + 1] + nc[i + 1, j] + nc[i, j + 1] + nc[i, j])
                w = rwgt[iwgt + 1]
                fgrad[1] += fx * w
                fgrad[2] += fy * w
            end
        end
        fgrad = abs.(fgrad)
        if fgrad[1] >= fgrad[2]
            jp, js = 1, 2
        else
            jp, js = 2, 1
        end
    end
    if check_dir >= 0 && check_dir + 1 != jp
        js, jp = jp, js
    end
    pdir[jp] = 1
    sdir[js] = 1

    if jp == 2
        f0[1, 2], f0[2, 1] = f0[2, 1], f0[1, 2]
    end

    if check_dir < 0
        vofi_check_secondary_side(impl_func, par, x0, h0, pdir, sdir, f0, xfs_pt, fth)
    end

    fx = (fc[3, 2] - fc[1, 2]) * have / (2 * hh[1])
    fy = (fc[2, 3] - fc[2, 1]) * have / (2 * hh[2])
    fxx = (fc[3, 2] + fc[1, 2] - 2 * fc[2, 2]) * have^2 / (hh[1]^2)
    fyy = (fc[2, 3] + fc[2, 1] - 2 * fc[2, 2]) * have^2 / (hh[2]^2)
    fxy = (fc[3, 3] - fc[3, 1] - fc[1, 3] + fc[1, 1]) * have^2 / (4 * hh[1] * hh[2])
    tmp = fx^2 + fy^2
    tmp = sqrt(tmp^3)
    Kappa = abs(fxx * fy^2 - 2 * fx * fy * fxy + fx^2 * fyy) / tmp
    a0, a1, a2, a3 = 2.30477, 28.5312, -46.2729, 56.9179
    est = a0 + Kappa * (a1 + Kappa * (a2 + a3 * Kappa))
    npt = Int(ceil(est))
    xfs_pt.ipt = max(4, min(npt, NGLM))
    return -1
end

function vofi_order_dirs_3D(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfsp)
    fc = @MArray zeros(vofi_real, NDIM, NDIM, NDIM)
    fd = @MMatrix zeros(vofi_real, NDIM, NDIM)
    hh = 0.5 .* h0
    fgrad = @MVector zeros(vofi_real, NDIM)
    np0 = 0
    nm0 = 0
    icc = -1
    check_dir = -1
    MIN_GRAD = 1.0e-4
    nmax0 = 8

    x = @MVector zeros(vofi_real, NDIM)
    for i in 0:1, j in 0:1, k in 0:1
        x[1] = x0[1] + i * h0[1]
        x[2] = x0[2] + j * h0[2]
        x[3] = x0[3] + k * h0[3]
        val = call_integrand(impl_func, par, x)
        f0[i + 1, j + 1, k + 1] = val
        if val > 0
            np0 += 1
        elseif val < 0
            nm0 += 1
        end
    end

    fgrad[1] = 0.25 * ((f0[2, 2, 2] + f0[2, 1, 2] + f0[2, 2, 1] + f0[2, 1, 1]) -
                       (f0[1, 2, 2] + f0[1, 1, 2] + f0[1, 2, 1] + f0[1, 1, 1])) / h0[1]
    fgrad[2] = 0.25 * ((f0[2, 2, 2] + f0[1, 2, 2] + f0[2, 2, 1] + f0[1, 2, 1]) -
                       (f0[2, 1, 2] + f0[1, 1, 2] + f0[2, 1, 1] + f0[1, 1, 1])) / h0[2]
    fgrad[3] = 0.25 * ((f0[2, 2, 2] + f0[2, 1, 2] + f0[1, 2, 2] + f0[1, 1, 2]) -
                       (f0[2, 2, 1] + f0[2, 1, 1] + f0[1, 2, 1] + f0[1, 1, 1])) / h0[3]
    fgradmod = max(sqrt(sum(abs2, fgrad)), MIN_GRAD)
    hm = maximum(hh)
    fth = sqrt(2.0) * fgradmod * hm

    if np0 * nm0 == 0
        n0 = @MArray zeros(Int, NSE, NSE, NSE)
        np0 = nm0 = 0
        for i in 0:1, j in 0:1, k in 0:1
            val = abs(f0[i + 1, j + 1, k + 1])
            if val > fth
                n0[i + 1, j + 1, k + 1] = 0
                if f0[i + 1, j + 1, k + 1] < 0
                    nm0 += 1
                else
                    np0 += 1
                end
            else
                n0[i + 1, j + 1, k + 1] = 1
            end
        end
        if nm0 == nmax0
            return 1
        elseif np0 == nmax0
            return 0
        end
        check_dir = vofi_check_boundary_surface(impl_func, par, x0, h0, f0, xfsp, n0)
        if check_dir < 0
            return nm0 > 0 ? 1 : 0
        end
    end

    for i in 0:1, j in 0:1, k in 0:1
        fc[2 * i + 1, 2 * j + 1, 2 * k + 1] = f0[i + 1, j + 1, k + 1]
    end
    x = @MVector zeros(vofi_real, NDIM)
    x[3] = x0[3] + hh[3]
    for i in 0:2:2, j in 0:2:2
        x[1] = x0[1] + i * hh[1]
        x[2] = x0[2] + j * hh[2]
        fc[i + 1, j + 1, 2] = call_integrand(impl_func, par, x)
    end
    x[1] = x0[1] + hh[1]
    x[2] = x0[2] + hh[2]
    for k in 0:2
        x[3] = x0[3] + k * hh[3]
        fc[2, 2, k + 1] = call_integrand(impl_func, par, x)
    end
    x[2] = x0[2] + hh[2]
    for i in 0:2:2
        x[1] = x0[1] + i * hh[1]
        for k in 0:2
            x[3] = x0[3] + k * hh[3]
            fc[i + 1, 2, k + 1] = call_integrand(impl_func, par, x)
        end
    end
    x[1] = x0[1] + hh[1]
    for j in 0:2:2
        x[2] = x0[2] + j * hh[2]
        for k in 0:2
            x[3] = x0[3] + k * hh[3]
            fc[2, j + 1, k + 1] = call_integrand(impl_func, par, x)
        end
    end

    fgrad .= 0
    for i in 0:1, j in 0:1, k in 0:1
        fx = 0.25 * ((fc[i + 2, j + 2, k + 2] + fc[i + 2, j + 1, k + 2] +
                      fc[i + 2, j + 2, k + 1] + fc[i + 2, j + 1, k + 1]) -
                     (fc[i + 1, j + 2, k + 2] + fc[i + 1, j + 1, k + 2] +
                      fc[i + 1, j + 2, k + 1] + fc[i + 1, j + 1, k + 1])) / hh[1]
        fy = 0.25 * ((fc[i + 2, j + 2, k + 2] + fc[i + 1, j + 2, k + 2] +
                      fc[i + 2, j + 2, k + 1] + fc[i + 1, j + 2, k + 1]) -
                     (fc[i + 2, j + 1, k + 2] + fc[i + 1, j + 1, k + 2] +
                      fc[i + 2, j + 1, k + 1] + fc[i + 1, j + 1, k + 1])) / hh[2]
        fz = 0.25 * ((fc[i + 2, j + 2, k + 2] + fc[i + 2, j + 1, k + 2] +
                      fc[i + 1, j + 2, k + 2] + fc[i + 1, j + 1, k + 2]) -
                     (fc[i + 2, j + 2, k + 1] + fc[i + 2, j + 1, k + 1] +
                      fc[i + 1, j + 2, k + 1] + fc[i + 1, j + 1, k + 1])) / hh[3]
        tmp = fc[i + 2, j + 2, k + 2] + fc[i + 2, j + 1, k + 2] +
              fc[i + 1, j + 2, k + 2] + fc[i + 1, j + 1, k + 2] +
              fc[i + 2, j + 2, k + 1] + fc[i + 2, j + 1, k + 1] +
              fc[i + 1, j + 2, k + 1] + fc[i + 1, j + 1, k + 1]
        tmp = max(abs(tmp), MIN_GRAD)
        fgrad[1] += fx / tmp
        fgrad[2] += fy / tmp
        fgrad[3] += fz / tmp
    end
    fgrad .= abs.(fgrad)
    jp = 1
    js = 2
    jt = 3
    if fgrad[1] >= fgrad[2]
        jp = 1
        js = 2
    else
        jp = 2
        js = 1
    end
    if fgrad[3] > fgrad[jp]
        jt = js
        js = jp
        jp = 3
    elseif fgrad[3] > fgrad[js]
        jt = js
        js = 3
    end

    if check_dir == 1 && xfsp[jp].isc[1] != 1
        if xfsp[js].isc[1] == 1
            jp, js = js, jp
        else
            jt, js, jp = js, jp, jt
        end
    end

    fill!(pdir, 0.0)
    fill!(sdir, 0.0)
    fill!(tdir, 0.0)
    pdir[jp] = 1.0
    sdir[js] = 1.0
    tdir[jt] = 1.0
    vofi_xyz2pst!(f0, jp, js, jt)

    if check_dir >= 0
        xfsp[5] = xfsp[jp]
    else
        vofi_check_secter_face(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfsp[5], fth)
    end
    vofi_check_tertiary_side(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfsp, fth)

    have = 0.5 * (h0[jp] + h0[js])
    curv = @MVector zeros(vofi_real, NDIM)
    sumf = @MVector zeros(vofi_real, NDIM)
    pd1, pd2, pd3 = Int(pdir[1]), Int(pdir[2]), Int(pdir[3])
    sd1, sd2, sd3 = Int(sdir[1]), Int(sdir[2]), Int(sdir[3])
    td1, td2, td3 = Int(tdir[1]), Int(tdir[2]), Int(tdir[3])
    for k in 1:NDIM
        kz = k - 1
        i0 = kz * td1
        j0 = kz * td2
        k0 = kz * td3
        sumf[k] = 0.0
        for i in 0:NDIM-1, j in 0:NDIM-1
            ii = i0 + i * sd1 + j * pd1
            jj = j0 + i * sd2 + j * pd2
            kk = k0 + i * sd3 + j * pd3
            fd[i + 1, j + 1] = fc[ii + 1, jj + 1, kk + 1]
            sumf[k] += fd[i + 1, j + 1]
        end
        sumf[k] = 1.0 / max(abs(sumf[k]), EPS_NOT0)
        fx = (fd[3, 2] - fd[1, 2]) * have / h0[js]
        fy = (fd[2, 3] - fd[2, 1]) * have / h0[jp]
        fxx = (fd[3, 2] + fd[1, 2] - 2 * fd[2, 2]) * have^2 / (hh[js] * hh[js])
        fyy = (fd[2, 3] + fd[2, 1] - 2 * fd[2, 2]) * have^2 / (hh[jp] * hh[jp])
        fxy = (fd[3, 3] - fd[3, 1] - fd[1, 3] + fd[1, 1]) * have^2 / (h0[js] * h0[jp])
        tmp = sqrt((fx^2 + fy^2)^3)
        curv[k] = abs(fxx * fy^2 - 2 * fx * fy * fxy + fx^2 * fyy) / tmp
    end
    Kappa = sum(sumf .* curv) / sum(sumf)
    a0, a1, a2, a3 = 2.34607, 16.5515, -5.53054, 54.0866
    tmp = a0 + Kappa * (a1 + Kappa * (a2 + a3 * Kappa))
    npt = Int(ceil(tmp))
    xfsp[5].ipt = max(4, min(npt, NGLM))
    return icc
end
function vofi_xyz2pst!(g0, jp, js, jt)
    if jt == 3 && js == 1
        g0[1, 2, 1], g0[2, 1, 1] = g0[2, 1, 1], g0[1, 2, 1]
        g0[1, 2, 2], g0[2, 1, 2] = g0[2, 1, 2], g0[1, 2, 2]
    elseif jp == 3
        g0[1, 1, 2], g0[2, 1, 1] = g0[2, 1, 1], g0[1, 1, 2]
        g0[1, 2, 2], g0[2, 2, 1] = g0[2, 2, 1], g0[1, 2, 2]
        if jt == 2
            g0[1, 1, 2], g0[1, 2, 1] = g0[1, 2, 1], g0[1, 1, 2]
            g0[2, 1, 2], g0[2, 2, 1] = g0[2, 2, 1], g0[2, 1, 2]
        end
    elseif js == 3
        g0[1, 1, 2], g0[1, 2, 1] = g0[1, 2, 1], g0[1, 1, 2]
        g0[2, 1, 2], g0[2, 2, 1] = g0[2, 2, 1], g0[2, 1, 2]
        if jp == 2
            g0[1, 1, 2], g0[2, 1, 1] = g0[2, 1, 1], g0[1, 1, 2]
            g0[1, 2, 2], g0[2, 2, 1] = g0[2, 2, 1], g0[1, 2, 2]
        end
    end
end
