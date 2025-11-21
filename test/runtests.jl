using Test
using VofiJul

@test vofi_real === Float64
@test NDIM == 3 == length(MinData().xval)

nodes3 = gauss_legendre_nodes(3)
weights4 = gauss_legendre_weights(4)
@test length(nodes3) == 3
@test nodes3[1] ≈ -0.7745966692414833
@test sum(weights4) ≈ 2.0

len = LenData()
@test length(len.xt0) == NGLM + 2

function plane_func(x, _)
    return x[1] + x[2] - 0.5
end

function neg_func(x, _)
    return -1.0
end

function pos_func(x, _)
    return 1.0
end

@test vofi_get_cell_type(neg_func, nothing, [0.0, 0.0], [1.0, 1.0], 2) == 1
@test vofi_get_cell_type(pos_func, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 3) == 0
@test vofi_get_cell_type(plane_func, nothing, [0.0, 0.0], [1.0, 1.0], 2) == -1

xex = zeros(Float64, 4)
cc_full = vofi_get_cc(neg_func, nothing, [0.0, 0.0], [1.0, 1.0], xex,
                      [0, 0], [0, 0], [0, 0], 2)
@test cc_full ≈ 1.0


function slanted_func(x, _)
    return 0.25 - x[1]
end

xex_cut2 = zeros(Float64, 4)
cc_cut2 = vofi_get_cc(slanted_func, nothing, [0.0, 0.0], [1.0, 1.0], xex_cut2,
                      [1, 1], [0, 0], [0, 0], 2)
@test cc_cut2 ≈ 0.75 atol=1e-6
@test xex_cut2[4] > 0.0

xex3 = zeros(Float64, 4)
cc_full3d = vofi_get_cc(neg_func, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], xex3,
                        [1, 0], [0, 0, 0, 0], [0, 0], 3)
@test cc_full3d ≈ 1.0
@test xex3[1:3] ≈ [0.5, 0.5, 0.5]

xex3_empty = zeros(Float64, 4)
cc_empty3d = vofi_get_cc(pos_func, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0],
                         xex3_empty, [0, 0], [0, 0, 0, 0], [0, 0], 3)
@test cc_empty3d ≈ 0.0

function plane_func3d(x, _)
    return x[1] + x[2] + x[3] - 1.5
end

xex3_plane = zeros(Float64, 4)
cc_plane = vofi_get_cc(plane_func3d, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0],
                       xex3_plane, [0, 0], [0, 0, 0, 0], [0, 0], 3)
@test cc_plane ≈ 0.5 atol=1e-3

# simple Cartesian integration of a sphere volume
function sphere_sdf(x, _)
    r = 0.4
    return sqrt((x[1])^2 + (x[2])^2 + (x[3])^2) - r
end

let
    n = 20
    h = 1.0 / n
    xmin = -0.5
    total_vol = 0.0
    cell_vol = h^3
    xex_tmp = zeros(Float64, 4)
    for i in 0:n-1, j in 0:n-1, k in 0:n-1
        xin = [xmin + i * h, xmin + j * h, xmin + k * h]
        cc = vofi_get_cc(sphere_sdf, nothing, xin, [h, h, h], xex_tmp,
                         [0, 0], [0, 0, 0, 0], [0, 0], 3)
        total_vol += cc * cell_vol
    end
    exact = 4 / 3 * π * 0.4^3
    @test total_vol ≈ exact atol=2e-2
end

# simple Cartesian integration of a circle area
function circle_sdf(x, _)
    r = 0.4
    return sqrt((x[1])^2 + (x[2])^2) - r
end

let
    n = 20
    h = 1.0 / n
    xmin = -0.5
    total_area = 0.0
    cell_area = h^2
    xex_tmp = zeros(Float64, 4)
    for i in 0:n-1, j in 0:n-1
        xin = [xmin + i * h, xmin + j * h]
        cc = vofi_get_cc(circle_sdf, nothing, xin, [h, h], xex_tmp,
                         [0, 0], [0, 0, 0, 0], [0, 0], 2)
        total_area += cc * cell_area
    end
    exact = π * 0.4^2
    @test total_area ≈ exact atol=2e-2
end