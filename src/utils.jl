# Utility functions for processing data and seismograms

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
    check_trace_order(t1, t2) -> a, b

Return the traces `t1` and `t2` as `a` and `b`, where `b` is always clockwise of `a`.
"""
check_trace_order(t1, t2) = Seis.angle_difference(t1.sta.azi, t2.sta.azi) > 0 ? (t1, t2) : (t2, t1)

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

"""
    maxabsloc2(A) -> imax, jmax

Return the indices of a two-dimensional array `A` where its maximum
absolute value occurs.
"""
function maxabsloc2(A)
    imin = jmin = 0
    maxval = zero(eltype(A))
    for j in 1:size(A, 2)
        for i in 1:size(A, 1)
            val = abs(A[i,j])
            if i == j == 1
                imin = jmin = 1
                maxval = val
            else
                if val > maxval
                    imin = i
                    jmin = j
                    maxval = val
                end
            end
        end
    end
    imin, jmin
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
