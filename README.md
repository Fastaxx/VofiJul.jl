# VofiJul.jl

Julia port of the [VOFI](https://github.com/vofi-dev/vofi) library for initializing volume fractions from an analytic implicit surface `f(x,y,z)=0`. Cells are treated as cuboids (2D or 3D), the reference phase is in the region `f < 0`, and integration follows the original VOFI algorithms.

## Usage (preview)

```julia
using VofiJul

# implicit function f(x) -> Float64 (negative inside)
sphere_sdf(x, _) = sqrt(sum(abs2, x)) - 0.4

xex = zeros(Float64, 4)                      # centroid / interface info
cc  = vofi_get_cc(sphere_sdf, nothing,
                  [-0.5, -0.5, -0.5],        # cell min corner
                  [1.0, 1.0, 1.0],           # cell sizes
                  xex, [1, 0], [0, 0, 0, 0], # centroids on
                  [0, 0], 3)                 # ndim=3
```

See `test/runtests.jl` for more examples, including a Cartesian integration of a sphere.

## Installation

The package is a pure Julia port; add it as a local dev package:

```shell
julia --project -e 'using Pkg; Pkg.develop(path=\".\"); Pkg.test()'
```

## Credits

All credits for the original algorithms, design, and C implementation go to the VOFI authors:
Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, Ruben Scardovelli, Philip Yecko, and Stéphane Zaleski. This Julia port is a reimplementation for the Julia ecosystem. 
