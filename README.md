# SeisSplit

SeisSplit is a package for measuring shear wave splitting using the
minimum eigenvalue method of Silver & Chan (1991), as modified by
Walsh et al. (2013).

It uses the [Seis.jl](https://github.com/anowacki/Seis.jl) package which
is an in-development community seismic analysis package in Julia.

## Installing

```julia
julia> ] # Press ']' to get to package mode

(v1.1)> add https://github.com/anowacki/Geodesics.jl https://github.com/anowacki/SAC.jl https://github.com/anowacki/Seis.jl https://github.com/anowacki/SeisSplit.jl
```

## Using

The `splitting` function performs shear wave splitting analysis and returns a named
tuple containing information about the analysis.  Provide two `Seis.Trace`s, and
optionaly specify the maximum delay time and number of fast orientaton and delay
time analysis points.

Example:

```julia
julia> using Seis, SeisSplit

julia> n, e = Seis.read_sac.(joinpath(dirname(pathof(SeisSplit)), "..", "test", "data", "wave.BH").*("N", "E"))
(Seis.Trace(.SWAV..BHN: delta=0.05, b=0.0, nsamples=1999), Seis.Trace(.SWAV..BHE: delta=0.05, b=0.0, nsamples=1999))

julia> s = splitting(n, e)
(phi = -90.0:1.0:90.0, dt = 0.1:0.1:4.0, lam1 = Float32[70.6024 71.7245 … 64.1079 64.2652; 70.5508 71.6332 … 64.6298 64.767; … ; 70.6525 71.813 … 63.5645 63.7426; 70.6024 71.7245 … 64.1079 64.2652], lam2 = Float32[7.24323 6.61557 … 10.8254 10.7367; 7.27189 6.66642 … 10.5352 10.4576; … ; 7.21534 6.56629 … 11.1276 11.0273; 7.24323 6.61557 … 10.8254 10.7367], phi_best = 40.0, dphi = 4.25, dt_best = 1.4, ddt = 0.1, spol = 10.013461f0, dspol = 1.1384997f0)

```

The following is the docstring for the `splitting` function:

```

    splitting(t1, t2; nphi=181, ndt=40, dt_max=4.0) -> results

Perform a search over a pair of Seis traces, `t1` and `t2`, for the smallest value of the
minimum eigenvalue of the covariance matrix between the traces, for a set of `nphi`×`ndt`
shear wave splitting operators, up to `dt_max` s.

`results` is a named tuple containing:

  •    `phi`: The set of fast shear wave orientations in °

  •    `dt`: The set of delays times in s

  •    `lam1`: The larger eigenvalues at each [phi,dt] point

  •    `lam2`: The smaller eigenvalues at each point

  •    `phi_best` and `dphi`: The best ϕ and the 1-sigma error therein

  •    `dt_best` and `ddt`: The best δt and the 1-sigma error therein

  •    `spol` and `dspol`: The source polarisation and an estimate of its uncertainty
      for the best-fitting ϕ and δt

```

## References

- Silver, P.G., Chan, W.W., 1991. Shear-wave splitting and subcontinental mantle
  deformation. J Geophys Res-Sol Ea 96, 16429–16454.
  doi:[10.1029/91JB00899](https://doi.org/10.1029/91JB00899)
- Walsh, E., Arnold, R., Savage, M.K., 2013. Silver and Chan revisited.
  Journal of Geophysical Research: Solid Earth 118, 5500–5515.
  doi[10.1002/jgrb.50386](https://doi.org/10.1002/jgrb.50386)
