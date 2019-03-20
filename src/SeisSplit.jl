"""
# SeisSplit

SeisSplit is a package for measuring shear wave splitting using the
minimum eigenvalue method of Silver & Chan (1991), as modified by
Walsh et al. (2013).

It uses the [Seis.jl](https://github.com/anowacki/Seis.jl) package which
is an in-development community seismic analysis package in Julia.


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

See the docstring for `splitting` for more information.


## References

- Silver, P.G., Chan, W.W., 1991. Shear-wave splitting and subcontinental mantle
  deformation. J Geophys Res-Sol Ea 96, 16429–16454.
  doi:[10.1029/91JB00899](https://doi.org/10.1029/91JB00899)
- Walsh, E., Arnold, R., Savage, M.K., 2013. Silver and Chan revisited.
  Journal of Geophysical Research: Solid Earth 118, 5500–5515.
  doi:[10.1002/jgrb.50386](https://doi.org/10.1002/jgrb.50386)
"""
module SeisSplit

import Distributions, FFTW
using LinearAlgebra
using StaticArrays

import Seis

export splitting


"Default number of δt points to search over"
const SPLIT_NDT = 40
"Default number of ϕ points to search over"
const SPLIT_NPHI = 181
"Default maximum δt (s)"
const SPLIT_DT_MAX = 4.0

"""
Struct containing the results of shear wave splitting analysis.
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
    "Origianl input trace 1"
    trace1::Seis.Trace{T,V}
    "Original input trace 2"
    trace2::Seis.Trace{T,V}
    "Analysis window start time (s)"
    window_start
    "Analysis windoe end time (s)"
    window_end
end

include("plots.jl")

"""
    splitting(t1, t2, window_start=starttime(t1), window_end=endtime(t1); nphi=$SPLIT_NPHI, ndt=$SPLIT_NDT, dt_max=$SPLIT_DT_MAX) -> results

Perform a search over a pair of Seis traces, `t1` and `t2`, for the smallest value of the
minimum eigenvalue of the covariance matrix between the traces, for a set of `nphi`×`ndt`
shear wave splitting operators, up to `dt_max` s.

`results` is a `SeisSplit.Result` containing:

- `phi`: The set of fast shear wave orientations in °
- `dt`: The set of delays times in s
- `lam1`: The larger eigenvalues at each [phi,dt] point
- `lam2`: The smaller eigenvalues at each point
- `phi_best` and `dphi`: The best ϕ and its 1σ uncertainty.  ϕ is measured
  clockwise from local north (towards east) in °.
- `dt_best` and `ddt`: The best δt and its 1σ uncertainty, in s
- `spol` and `dspol`: The source polarisation and an estimate of its uncertainty for the
  best-fitting ϕ and δt.  `spol` is given in ° clockwise of local north.
- `trace1` and `trace2`, the original input traces, where `trace2` is clockwise of `trace1`
- `window_start`, `window_end`, the analysis time window end points.
"""
function splitting(t1::Seis.Trace{T,V}, t2::Seis.Trace{T,V},
                   window_start=Seis.starttime(t1), window_end=Seis.endtime(t1);
                   nphi=SPLIT_NPHI, ndt=SPLIT_NDT, dt_max=SPLIT_DT_MAX) where {T,V}
    t1, t2 = check_trace_order(t1, t2) # t2 is now clockwise of t1
    n, e = Seis.rotate_through(t1, t2, -t1.sta.azi)
    Seis.cut!.((n, e), window_start, window_end)
    phi = range(-90., stop=90., length=nphi)
    dt = range(dt_max/ndt, stop=dt_max, length=ndt)
    lam1, lam2 = zeros(T, nphi, ndt), zeros(T, nphi, ndt)
    N, E = deepcopy(n), deepcopy(e)
    @inbounds for j in eachindex(dt), i in eachindex(phi)
        lam1[i,j], lam2[i,j] = compute_eigvals(n, e, phi[i], dt[j], N, E)
    end
    # Best parameters
    ip, idt = minloc2(lam2)
    phi_best = phi[ip]
    dt_best = dt[idt]
    # Find source polarisation and error therein
    spol, dspol = sourcepol(n, e, phi[ip], dt[idt])
    # Errors in splitting parameters
    dphi, ddt = errors_scale_lam2!(lam2, N, E, n, e, phi_best, dt_best, spol, phi, dt)
    Result(phi, dt, lam1, lam2, phi_best, dphi, dt_best, ddt, spol, dspol, t1, t2,
           window_start, window_end)
end

"""
    check_trace_order(t1, t2) -> a, b

Return the traces `t1` and `t2` as `a` and `b`, where `b` is always clockwise of `a`.
"""
check_trace_order(t1, t2) = Seis.angle_difference(t1.sta.azi, t2.sta.azi) > 0 ? (t1, t2) : (t2, t1)

"""
    apply_split!(n, e, phi, dt)

Apply a splitting operator to two traces.  `n` is assumed to be oriented north,
and `e` east.
"""
function apply_split!(n::T, e::T, phi::Number, dt::Number) where {T<:Seis.Trace}
    rotate_traces!(n, e, -phi)
    # Move the fast trace forward
    shift_trace!(n, -dt)
    rotate_traces!(n, e, phi)
    nothing
end

"""
    rotate_traces!(s1, s2, phi)

Rotate the pair of traces `s1` and `s2` clockwise through `phi`°.  Particle motion
appears to move anticlockwise.
"""
function rotate_traces!(s1, s2, phi)
    # Minimal version of Seis.rotate_through! which does not update headers
    sinp, cosp = sincos(deg2rad(phi))
    @inbounds for i in 1:Seis.nsamples(s1)
        s1.t[i], s2.t[i] = cosp*s1.t[i] - sinp*s2.t[i], sinp*s1.t[i] + cosp*s2.t[i]
    end
end

"""
    shift_trace!(t, dt) -> t

Shift a trace `t` backwards in time by `dt` seconds.
"""
function shift_trace!(t, dt)
    T = eltype(t.t)
    n = round(Int, dt/t.delta)
    if n == 0
        return t
    end
    array_shift!(t.t, n)
end

"""
    array_shift!(a, n, pad=zero(eltype(a))) -> a

Shift the values within an array back by `n` places, and replace
the shifted values with `pad`.  This function
works similarly to `circshift`, but rather
than elements shifting back round to the other end of the array when
shifted off one end, replace the shifted elements with `pad`.
"""
function array_shift!(a, n, pad=zero(eltype(a)))
    n == 0 && return a
    N = length(a)
    # Just zero out trace if shifted the whole way
    if abs(n) >= N
        a .= pad
        return a
    end
    if n > 0
        @inbounds for i in N:-1:(n+1)
            a[i] = a[i-n]
        end
        @inbounds a[1:n] .= pad
    elseif n < 0
        @inbounds for i in 1:(N+n)
            a[i] = a[i-n]
        end
        @inbounds a[end+n+1:end] .= pad
    end
    a
end

"""
    compute_eigvals(t1, t2, phi, dt, T1, T2) -> λ₁, λ₂

Compute the two eigenvalues for traces for a specific fast orientation, `phi`°,
and delay time, `dt` s.  `t1` and `t2` are assumed to be orthogonal, and `T1`
and `T2` are a pair of temporary storage traces matching `t1` and `t2`, respectively.
"""
function compute_eigvals(t1::Seis.Trace, t2::Seis.Trace, phi::Number, dt::Number,
                         T1::Seis.Trace, T2::Seis.Trace)
    # Note that we swap the delay time, because we are applying the inverse
    # splitting operator.
    copy_trace!(T1, t1)
    copy_trace!(T2, t2)
    apply_split!(T1, T2, phi, -dt)
    c = covar(T1, T2)
    eval = eigvals(c)
    maximum(eval), minimum(eval)
end

"""
    sourcepol(t1, t2, phi, dt) -> spol, dspol

Compute the source polarisation `spol` and error `dspol` for the pair of traces `t1`
and `t2` for the best-fitting splitting parameters `phi`° and `dt` s.
"""
function sourcepol(t1, t2, phi::Number, dt::Number)
    T1, T2 = deepcopy(t1), deepcopy(t2)
    apply_split!(T1, T2, phi, -dt)
    c = covar(T1, T2)
    eval, evec = eigen(c)
    i1 = argmax(eval)
    i2 = 3 - i1
    lam1 = eval[i1]
    lam2 = eval[i2]
    spol = mod(rad2deg(atan(evec[2,i1], evec[1,i1])), 180)
    dspol = rad2deg(atan(lam2/lam1))
    spol, dspol
end

"""
    covar(t1, t2) -> c

Return the covariance matrix `c` for two traces.
"""
function covar(t1::Seis.Trace{T,V,S}, t2::Seis.Trace{T,V,S}) where {T,V,S}
    c11 = c22 = c12 = zero(eltype(V))
    @inbounds for i in 1:Seis.nsamples(t1)
        c11 += t1.t[i]^2
        c22 += t2.t[i]^2
        c12 += t1.t[i]*t2.t[i]
    end
    c = SArray{Tuple{2,2},T}(c11, c12, c12, c22)
    c
end

"""
    errors_scale_lam2!(lam2, N, E, n, e, phi_best, dt_best, spol, phi, dt) -> σϕ, σδt

Return the 1σ errors in ϕ and δt, `σϕ` and `σδt`, for the best-fitting splitting parameters
`phi_best` and `dt_best`, and scale the λ₂ surface, `lam2`, such that the 95% confidence contour
has a value of 1.  `spol` is the retrieved source polarisation.
`n` and `e` are the north and east input traces, and `N` and `E` are
temporary traces which will be rotated to the fast-slow direction to compute
the errors.  `phi` and `dt` are the values of ϕ and δt along the axes of `lam2`.
"""
function errors_scale_lam2!(lam2, N, E, n, e, phi_best, dt_best, spol, phi, dt)
    copy_trace!(N, n)
    copy_trace!(E, e)
    # Rotate to fast-slow and remove splitting
    rotate_traces!(N, E, -phi_best)
    shift_trace!(N, dt_best)
    # Rotate to polarisation direction.  E should be ~ 0 and thus represent noise
    rotate_traces!(N, E, phi_best - spol)
    ν = degrees_of_freedom(E.t)
    # Critical F-statistic at 95% F(k, ν-2) where k = 2
    Fν = Distributions.quantile(Distributions.FDist(2, ν - 2), 0.95)
    # Scale λ₂
    λ₂_min = minimum(lam2)
    if λ₂_min < eps(λ₂_min)
        @warn("Minimum λ₂ is zero (can occur in synthetics).  " *
              "Normalising by non-zero value.  σ may be inaccurate.")
        λ₂_min = √eps(λ₂_min)
    end
    # Normalise 95% contour to 1
    lam2 .= lam2./(1 + 2Fν/(ν - 2))./λ₂_min
    # Get errors by scanning the error surface
    σ_ϕ, σ_δt = errors_from_normalised_lam2(lam2, phi, dt)
end

"""
    degrees_of_freedom(f) -> ndf

Return an estimation of the degrees of freedom in the noise vector `f`.
"""
function degrees_of_freedom(f)
    # Compute amplitude spectrum
    F = abs.(FFTW.rfft(f))
    E₂ = sum(x->x^2, F)
    E₂ -= (F[1]^2 + F[end]^2)/2
    E₄ = 1/3*(F[1]^4 + F[end]^4)
    for i in 2:(length(F)-1)
        E₄ += 4/3*F[i]^4
    end
    ndf = round(Int, 4E₂^2/E₄ - 2)
end

"""
    errors_from_normalised_lam2(lam2, phi, dt) -> σ_ϕ, σ_δt

Compute the 1σ uncertainties in ϕ and δt for the normalised λ₂ surface `lam2`,
given the values of ϕ `phi` and δt `dt` along the axes of `lam2`.
"""
function errors_from_normalised_lam2(lam2, phi, dt)
    I = findall(x->x<=1, lam2)
    idt1 = minimum(x->x[2], I)
    idt2 = maximum(x->x[2], I)
    σ_δt = (idt2 - idt1)*step(dt)/4
    nphi = length(unique(getindex.(I, 1)))
    σ_ϕ = nphi*step(phi)/4
    σ_ϕ, σ_δt
end

"""
    copy_trace!(to::T, from::T) where {T<:Seis.Trace}

Copy the `Seis.Trace` trace `a` into `b` without allocating.
"""
@inline function copy_trace!(to::T, from::T) where {T<:Seis.Trace}
    to.t .= from.t
    nothing
end

"""
    minloc2(a) -> imin, jmin

Return an 2-tuple which contains the location `(imin,jmin)` of the minimum value of a
2-dimensional array.
"""
function minloc2(A)
    imin = jmin = 0
    min = zero(eltype(A))
    for j in 1:size(A, 2)
        for i in 1:size(A, 1)
            if i == j == 1
                imin = jmin = 1
                min = A[i,j]
            else
                if A[i,j] < min
                    imin = i
                    jmin = j
                    min = A[i,j]
                end
            end
        end
    end
    imin, jmin
end

end # module
