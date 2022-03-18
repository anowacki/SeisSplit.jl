"""
# SeisSplit

SeisSplit is a package for measuring shear wave splitting using the
minimum eigenvalue method of Silver & Chan (1991), as modified by
Walsh et al. (2013).

It uses the [Seis.jl](https://github.com/anowacki/Seis.jl) package which
is an in-development community seismic analysis package in Julia.

Computation of splitting using the rotation-correlation method
(e.g., Bowman & Ando, 1987) is also possible.  This is computed
automatically by the main `splitting` function for quality control
purposes.

## Using

The `splitting` function performs shear wave splitting analysis and returns a `SeisSplit.Result`
type containing information about the analysis.  Provide two `Seis.Trace`s, and
optionally specify the maximum delay time and number of fast orientation and delay
time analysis points.

Example:

```julia
julia> using Seis, SeisSplit

julia> n, e = Seis.read_sac.(joinpath(dirname(pathof(SeisSplit)), "..", "test", "data", "wave.BH").*("N", "E"))
(Seis.Trace(.SWAV..BHN: delta=0.05, b=0.0, nsamples=1999), Seis.Trace(.SWAV..BHE: delta=0.05, b=0.0, nsamples=1999))

julia> s = splitting(n, e)
SeisSplit.Result{Float32,Array{Float32,1}}(-90.0:1.0:90.0, 0.0:0.10256410256410256:4.0, Float32[69.407 70.6024 … 64.1079 64.2652; 69.407 70.5508 … 64.6298 64.767; … ; 69.407 70.6525 … 63.5645 63.7426; 69.407 70.6024 … 64.1079 64.2652], Float32[8.79587 8.05464 … 12.038 11.9395; 8.79587 8.08652 … 11.7154 11.6291; … ; 8.79587 8.02363 … 12.3742 12.2627; 8.79587 8.05464 … 12.038 11.9395], 41.0, 0.5, 1.3333333333333333, 0.0, 10.455532f0, 1.1246407f0, Seis.Trace(.SWAV..BHN: delta=0.05, b=0.0, nsamples=1999), Seis.Trace(.SWAV..BHE: delta=0.05, b=0.0, nsamples=1999), 0.0f0, 99.9f0)

```
See the docstring for `splitting` for more information.


### Plotting results

You can create a diagnostic plot of a `SeisSplit.Result` by loading
[`Plots.jl`](https://github.com/JuliaPlots/Plots.jl) and calling `plot()` on the result:

```julia
julia> using Plots

julia> plot(s)
```


## References

- Bowman, J. R., Ando, M., 1987. Shear-wave splitting in the
  upper-mantle wedge above the Tonga subduction zone.  Geophys J R
  Astr Soc 88, 25–41.
  doi:[10.1111/j.1365-246X.1987.tb01367.x](https://doi.org/10.1111/j.1365-246X.1987.tb01367.x)
- Silver, P.G., Chan, W.W., 1991. Shear-wave splitting and subcontinental mantle
  deformation. J Geophys Res-Sol Ea 96, 16429–16454.
  doi:[10.1029/91JB00899](https://doi.org/10.1029/91JB00899)
- Walsh, E., Arnold, R., Savage, M.K., 2013. Silver and Chan revisited.
  Journal of Geophysical Research: Solid Earth 118, 5500–5515.
  doi:[10.1002/jgrb.50386](https://doi.org/10.1002/jgrb.50386)

"""
module SeisSplit

import DelimitedFiles, DSP, Printf, Statistics

import Distributions, FFTW
using LinearAlgebra
using StaticArrays

import Seis

export
    lam_ratio,
    lam1_ratio,
    lam2_ratio,
    quality,
    rotation_correlation,
    snr_restivo_helffrich,
    splitting


"Default number of δt points to search over"
const SPLIT_NDT = 41
"Default number of ϕ points to search over"
const SPLIT_NPHI = 181
"Default maximum δt (s)"
const SPLIT_DT_MAX = 4.0

"""
  Result{T,V}

Struct containing the results of shear wave splitting analysis.

Note that the two trace fields, `trace1` and `trace2` may not be in the same order
as originally given to the [`splitting`](@ref) function.  This is because
they are re-ordered so that `trace2` is clockwise of `trace1`.  Hence the
frame of reference in which this result is measured is determined by the
two traces' orientations.
"""
struct Result{T,V}
    "Range of values of ϕ over which we searched (° clockwise from north)"
    phi
    "Values of δt searched (s)"
    dt
    "Larger eigenvalue at each point in (ϕ,δt) space (normalised)"
    lam1
    "Smaller eigenvalue"
    lam2
    "Optimal fast orientation ϕ (° clockwise from north)"
    phi_best
    "Uncertainty in optimal ϕ (°)"
    dphi
    "Optimal δt (s)"
    dt_best
    "Uncertainty in optimal δt (s)"
    ddt
    "Source polarisation for optimal splitting (° clockwise from north)"
    spol
    "Uncertainty in source polarisation"
    dspol
    "Original input trace 1"
    trace1::Seis.Trace{T,V}
    "Original input trace 2"
    trace2::Seis.Trace{T,V}
    "Analysis window start time (s)"
    window_start
    "Analysis window end time (s)"
    window_end
    "Degrees of freedom"
    ndf
    "Cross correlation range of ϕ values (° clockwise from north)"
    xcorr_phi
    "Cross correlation range of δt values searched (s)"
    xcorr_dt
    "Cross correlation at each point in (ϕ,δt) space"
    xcorr_map
    "Maximum cross correlation point in ϕ (°)"
    xcorr_phi_best
    "Maximum cross correlation point in δt (s)"
    xcorr_dt_best
    "Frame of reference for ϕ: `:geographic` for azimuth from north; `:trace` for
    `trace1` to `trace2`"
    reference_frame
end

include("utils.jl")
include("io.jl")
include("qc.jl")
include("minimum_eigenvalue.jl")
include("rotation_correlation.jl")
include("plots.jl")

end # module
