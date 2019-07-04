# Functions for the rotation-correlation method

"""
    rotation_correlation(t1, t2, window_start=starttime(t1), window_end=endtime(t1); nphi=$SPLIT_NPHI, dt_max=$SPLIT_DT_MAX) -> dt, phi, xcorr_map, phi_best, dt_best

Compute the rotation-correlation map between two traces `t1` and `t2`

`phi` is a set of orientations to test, whilst `dt_max` is the largest delay
time tested.

`xcorr_map[idt,iphi]` holds the normalised cross correlation at `dt[idt]`
 and `phi[iphi]`.  `phi_best` and `dt_best` give the value of the fast
 orientation and delay time, respectively, at the maximum correlation point.
"""
function rotation_correlation(t1::Seis.Trace{T,V}, t2::Seis.Trace{T,V},
                              window_start=Seis.starttime(t1), window_end=Seis.endtime(t1);
                              nphi=SPLIT_NPHI, dt_max=SPLIT_DT_MAX) where {T,V}
    t1, t2 = check_trace_order(t1, t2) # t2 is now clockwise of t1
    n, e = Seis.rotate_through(t1, t2, -t1.sta.azi)
    Seis.cut!.((n, e), window_start, window_end)
    phi = range(-90., stop=90., length=nphi+1)[1:nphi]
    _rotation_correlation!(n, e, phi, dt_max)
end

"""
    _rotation_correlation!(n, e, phi, dt_max) -> dt, phi, xcorr_map, phi_best, dt_best

Compute the rotation-correlation map between the two traces `n` and `e`.
`n` must point northwards and `e` eastwards.

`phi` is a set of orientations to test, whilst `dt_max` is the largest delay
time tested.

`xcorr_map[idt,iphi]` holds the normalised cross correlation at `dt[idt]`
 and `phi[iphi]`.  `phi_best` and `dt_best` give the value of the fast
 orientation and delay time, respectively, at the maximum correlation point.
"""
function _rotation_correlation!(n::Seis.Trace{T,V}, e::Seis.Trace{T,V}, phi, dt_max) where {T,V}
    npts = Seis.nsamples(n)
    Seis.nsamples(e) == npts || throw(ArgumentError("traces must be the same length"))
    nphi = length(phi)
    nxcorr = 2npts - 1
    # Correlation map in dt and phi
    xcorr_map = Array{T}(undef, nxcorr, nphi)
    xcorr = Vector{T}(undef, nxcorr)
    # Times for each cross-correlation point
    time = ((1-npts):(npts-1)).*n.delta
    # Create work arrays and FFT plan for (I)FFTs.  This is just a tuple of stuff to
    # splat into the call to _xcorr_fd!()
    work_arrays = _correlation_arrays(T, npts)
    # Calculate correlation map
    ϕ_prev = 0.0
    @inbounds for (iϕ, ϕ) in enumerate(phi)
        δϕ = ϕ - ϕ_prev
        rotate_traces!(n, e, -δϕ)
        _xcorr_fd!(xcorr, Seis.trace(n), Seis.trace(e), work_arrays...; norm=true)
        xcorr_map[:,iϕ] .= xcorr
        ϕ_prev = ϕ
    end
    # Orient map to correct location
    idt, iphi = maxabsloc2(xcorr_map)
    dt_best = abs(time[idt])
    # Flip map if fast and slow anticorrelated
    xcorr_map[idt,iphi] < 0 && (xcorr_map .*= -1)
    xcorr_map = if time[idt] >= 0
        # Fast wave arrives later, so take positive time half
        # of map and then rotate map 90°
        iphi = mod1(iphi + nphi÷2, nphi)
        xcorr_map = xcorr_map[0 .<= time .<= dt_max, :]
        circshift(xcorr_map, (0, nphi÷2))
    else
        # Fast wave arrives first, so take negative time half
        # of map
        xcorr_map = reverse(xcorr_map[-dt_max .<= time .<= 0, :], dims=1)
    end
    phi_best = phi[iphi]
    time[0 .<= time .<= dt_max], phi, xcorr_map, phi_best, dt_best
end

"""
    _correlation_arrays(T, n) -> (apad, bpad, Apad, Bpad, conv, fplan, bplan)

Return a tuple of work arrays to be used be `_xcorr_fd!` to perform
in-place cross-correlation of two traces of length `n` and element
type `T`.
"""
function _correlation_arrays(T, n)
    nxc = 2n - 1
    np2 = nxc > 1024 ? nextprod([2,3,5], nxc) : nextpow(2, nxc)
    np = np2÷2 + 1
    # Input vectors padded to power-of-2 length
    apad = Vector{T}(undef, np2)
    bpad = similar(apad)
    # Fourier transformed vectors
    Apad = Vector{Complex{T}}(undef, np)
    Bpad = similar(Apad)
    # Convolution of two vectors
    conv = similar(Apad)
    # FFTW transform plans
    fplan = FFTW.plan_rfft(apad)
    bplan = FFTW.plan_irfft(conv, np2)
    (apad, bpad, Apad, Bpad, conv, fplan, bplan)
end

"""
    _xcorr_fd!(xc, a, b, apad, bpad, Apad, Bpad, conv, fplan, bplan; norm=false) -> xc

Compute the cross correlation between arrays `a` and `b`, which must
have the same length.  `xc` is filled in-place and returned.

Remaining arguments are work arrays and FFT forward and backward
real-to-complex plans computed with a call to `_correlation_arrays`.

If `norm` is `true`, then the cross correlation is normalised by
`sqrt(ACF(a) * ACF(b))`, where `ACF` is the autocorrelation.

For speed, no bounds checks are performed in this function and
therefore behaviour is undefined if any arrays are unexpected sizes.
"""
function _xcorr_fd!(xc, a, b, apad, bpad, Apad, Bpad, conv, fplan, bplan; norm=false)
    n = length(a)
    nxc = 2n - 1
    @inbounds apad[1:n] .= a
    @inbounds for i in 1:n
        bpad[i] = b[end-i+1]
    end
    @inbounds apad[(n+1):end] .= 0
    @inbounds bpad[(n+1):end] .= 0
    LinearAlgebra.mul!(Apad, fplan, apad)
    LinearAlgebra.mul!(Bpad, fplan, bpad)
    @inbounds conv .= Apad .* Bpad
    LinearAlgebra.mul!(apad, bplan, conv)
    @inbounds xc .= apad[1:nxc]
    if norm
        norm_val = sum(x->x^2, a) * sum(x->x^2, b)
        norm_val = √norm_val
        xc ./= norm_val
    end
    xc
end

"""
    cross_correlation_td!(xc, a, b) -> xc

Fill the contents of `xc` with the cross-correlation of the two vectors `a` and `b`,
which must be the same length `n`.  `xc` must have length `2n - 1`.  The zero-lag
point is therefore at `xc[n]`.

This is done in the time domain; rough benchmarking suggests that this can be faster
than a frequency-domain approach when `a` and `b` are less than ~350 in length.
"""
function cross_correlation_td!(xc, a, b)
    n = length(a)
    @boundscheck length(b) == n ||
        throw(ArgumentError("two vectors must be the same length"))
    @boundscheck length(xc) == 2n - 1 ||
        throw(ArgumentError("output vector must be `2n - 1` long"))
    xc .= 0
    i = 0
    @inbounds for ilag in (1-n):(n-1)
        i += 1
        iover = n - abs(ilag)
        if ilag <= 0
            @simd for j in 1:iover
                xc[i] += a[j]*b[j-ilag]
            end
        else
            @simd for j in 1:iover
                xc[i] += a[j+ilag]*b[j]
            end
        end
    end
    xc
end

"""
    cross_correlation_td(a, b) -> xc

Return the cross-correlation of the two vectors `a` and `b`, which must be the same
length.

This is done in the time domain; rough benchmarking suggests that this can be faster
than a frequency-domain approach when `a` and `b` are less than ~350 in length.
"""
cross_correlation_td(a, b) = cross_correlation!(similar(a, 2length(a)-1), a, b)

"""
    cross_correlation_fd!(xc, a, b) -> xc

Fill `xc` with the cross-correlation of the two vectors `a` and `b`, which must be
the same length `n`.  `xc` must have length `2n - 1`.  This is done in the frequency domain.
"""
cross_correlation_fd!(xc, a, b) = (xc .= DSP.xcorr(a, b); xc)

"""
    cross_correlation_fd(a, b) -> xc

Return the cross-correlation of the two vectors `a` and `b`, which must be the same
length.  This is done in the frequency domain.
"""
cross_correlation_fd(a, b) = DSP.xcorr(a, b)

"""
    cross_correlation!(xc, a, b, cutoff=350) -> xc

Fill `xc` with the cross-correlation between the vectors `a` and `b`, which must be
the same length `n`.  `xc` must have length `2n - 1`.

When `n` is less than `cutoff`, a time-domain implementation is used, whilst
when `n` is `cutoff` or larger, a frequency-domain implementation is used.
"""
cross_correlation!(xc, a, b, cutoff=350) = length(a) < cutoff ?
                                               cross_correlation_td!(xc, a, b) :
                                               cross_correlation_fd!(xc, a, b)

"""
    cross_correlation(a, b, cutoff=350) -> xc

Return the cross-correlation  `xc` between the vectors `a` and `b`, which must be the same
length `n` and have the same element type.  `xc` has length `2n - 1`.

When `n` is less than `cutoff`, a time-domain implementation is used, whilst
when `n` is `cutoff` or larger, a frequency-domain implementation is used.
"""
cross_correlation(a::AbstractVector{T}, b::AbstractVector{T}, cutoff=350) where T =
    length(a) < cutoff ?
        cross_correlation_td!(Vector{T}(undef, 2*length(a) - 1), a, b) :
        cross_correlation_fd(a, b)
