module SeisSplitMakieExt

import DSP
import Makie
import Seis
using SeisSplit: SeisSplit, Result

function __init__()
    "Loaded SeisSplitMakieExt"
end

SeisSplit.plot_result(s::Result; kwargs...) = Makie.plot(s; kwargs...)

"""
    plot!(fig, s::SeisSplit.Result; max_samples=1000, antialias=true) -> ::Plots.Plot

Plot the results of shear wave splitting analysis from the `SeisSplit.splitting` function.

# Keyword arguments
- `antialias = true`: Control whether decimation of traces for plotting purposes
  uses antialising to avoid spurious signals.  `true` will avoid spurious aliasing
  signals, whilst `false` may be quicker.
- `max_samples = 1000`: Maximum number of samples to plot in the trace and particle
  motion subplots.  Windows or traces with more than this number are decimated to
  below this number to speed up plotting.
"""
function Makie.plot(s::Result; max_samples=1000, antialias=true, axis=())
    # Common axis options
    axis = (; xgridvisible=false, ygridvisible=false, axis...)

    # Input traces and whether we are now in the N-E frame
    n_orig, e_orig, traces_n_e = if all(Seis.is_horizontal, (s.trace1, s.trace2))
        # Rotated to N-E if horizontal
        (Seis.rotate_through(s.trace1, s.trace2, -s.trace1.sta.azi)..., true)
    else
        # Otherwise just in the input orientation and hence ϕ is from
        # trace1 to trace2
        (s.trace1, s.trace2, false)
    end

    # Tuple of window limits
    window = s.window_start, s.window_end

    fig = Makie.Figure(; size=(800, 700))

    ## Input traces in original orientation with analysis window highlighted
    # Traces
    trace1, trace2 = _decimate_to_max_samples.((s.trace1, s.trace2), max_samples,
        antialias)
    ax1, lp1 = Makie.lines(
        fig[1,1],
        Seis.times(trace1),
        Seis.trace(trace1);
        label=s.trace1.sta.cha,
        axis=(title="Input traces", axis...),
    )
    Makie.lines!(ax1, Seis.times(trace2), Seis.trace(trace2); label=s.trace2.sta.cha)
    Makie.axislegend(ax1)

    # Window limits
    Makie.vlines!(ax1, [window...]; color=:red, label="")

    ## Rotated to spol
    spol, spol_90 = _decimate_to_max_samples.(
        Seis.rotate_through(n_orig, e_orig, #=180-s.spol=#s.spol),
        max_samples, antialias, Ref(window)
    )

    ax2, lp2 = Makie.lines(
        fig[1,2], Seis.times(spol), Seis.trace(spol);
        color=:blue,
        axis=(
            title="Source pol. pre-corr",
            limits=(extrema(Seis.times(spol)), nothing),
            axis...
        )
    )
    Makie.lines!(ax2, Seis.times(spol_90), Seis.trace(spol_90); color=:red)

    ## Corrected and rotated to spol
    n, e = _decimate_to_max_samples.((n_orig, e_orig), max_samples, antialias, Ref(window))
    SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
    fast, slow = Seis.rotate_through(n, e, #=180-s.spol=#s.spol)
    ax3, pl3 = Makie.lines(
        fig[1,3], Seis.times(fast), Seis.trace(fast);
        color=:blue,
        axis=(
            title="Source pol. corr",
            limits=(extrema(Seis.times(fast)), nothing),
            axis...
        )
    )
    Makie.lines!(ax3, Seis.times(slow), Seis.trace(slow); color=:red, label="Slow")

    ## Fast and slow before dt shift
    sfast, sslow = Seis.rotate_through!.(_decimate_to_max_samples.((n_orig, e_orig), max_samples,
        antialias, Ref(window))..., s.phi_best)
    # Flip slow wave if opposite polarity to make it easier to see if traces line up
    slow_polarity = let p = _apparent_polarity(sfast, sslow)
        # Account for zero polarity if traces not correlated at all
        p == 0 ? 1 : p
    end
    ax4, pl4 = Makie.lines(
        fig[2,1], Seis.times(sfast), Seis.trace(sfast);
        color=:blue, label="Fast",
        axis=(limits=(window, nothing), title="Fast-slow pre-corr", axis...),
    )
    Makie.lines!(ax4, Seis.times(sfast), slow_polarity.*Seis.trace(sslow); color=:red, label="Slow")

    ## Fast and slow after dt shift
    n, e = _decimate_to_max_samples.((n_orig, e_orig), max_samples, antialias, Ref(window))
    # Note we are taking the splitting off, hence -δt
    SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
    fast, slow = Seis.rotate_through(n, e, s.phi_best)
    ax5, pl5 = Makie.lines(
        fig[2,2], Seis.times(fast), Seis.trace(fast);
        color=:blue,
        axis=(limits=(window, nothing), title="Fast-slow corr", axis...)
    )
    Makie.lines!(ax5, Seis.times(fast), slow_polarity.*Seis.trace(slow); color=:red)

    ## λ₂ surface
    # 10σ - 100σ contours
    ax6, pl6 = Makie.contour(
        fig[2,3], s.dt, s.phi, s.lam2';
        levels=10:10:100, linewidth=1, color=(:gray, 0.6),
        axis=(
            xlabel="δt / s",
            ylabel="ϕ / °",
            limits=(extrema(s.dt), extrema(s.phi)),
            title="λ₂",
            axis...
        ),
    )
    # 1σ – 10σ contours
    Makie.contour!(
        ax6, s.dt, s.phi, s.lam2';
        levels=1:10, linewidth=1, color=:black
    )
    # 1σ contour
    Makie.contour!(
        ax6, s.dt, s.phi, s.lam2';
        levels=[1], linewidth=3, color=:black
    )
    # Best value
    Makie.scatter!(ax6, s.dt_best, s.phi_best)

    ## Particle motion before correction
    amax = 1.2*maximum(x -> sqrt(x[1]^2 + x[2]^2), zip(Seis.trace(e_orig), Seis.trace(n_orig)))
    x, y = _decimate_to_max_samples.(
        (e_orig, n_orig), max_samples, antialias, Ref(window)
    )
    ax7, pl7 = Makie.lines(
        fig[3,1], Seis.trace(x), Seis.trace(y);
        color=:black,
        axis=(
            aspect=Makie.DataAspect(),
            limits=(-amax, amax, -amax, amax),
            xlabel=(traces_n_e ? "East" : e_orig.sta.cha),
            ylabel=(traces_n_e ? "North" : n_orig.sta.cha),
            title="PM pre-corr",
            axis...
        ),
    )

    ## Particle motion after correction
    n, e = deepcopy.((n_orig, e_orig))
    SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
    y, x = _decimate_to_max_samples.((n, e), max_samples, antialias, Ref(window))

    ax8, pl8 = Makie.lines(
        fig[3,2], Seis.trace(x), Seis.trace(y);
        color=:black,
        axis=(
            aspect=Makie.DataAspect(),
            limits=(-amax, amax, -amax, amax),
            xlabel=(traces_n_e ? "East" : e_orig.sta.cha),
            ylabel=(traces_n_e ? "North" : n_orig.sta.cha),
            title="PM corr",
            axis...
        ),
    )

    fig
end

"""
    _decimate_to_max_samples(t, max_samples, antialias[, cut]) -> t′

Return a new trace `t′` cut between `cut[1]` and `cut[2]`, and decimated such that
its number of samples is less than or equal to `max_samples`.  If `antialias`
is `true`, then an antialiasing filter is applied; otherwise it is not.
"""
function _decimate_to_max_samples(t, max_samples, antialias, cut=(Seis.starttime(t), Seis.endtime(t)))
    max_samples > 10 || throw(ArgumentError("`max_samples` must be 10 or more"))
    t_cut = Seis.cut(t, cut...)
    n = Seis.nsamples(t_cut)
    decimation = 1
    while n/decimation > max_samples
        decimation += 1
    end
    Seis.decimate!(t_cut, decimation, antialias=antialias)
end

"""
    _apparent_polarity(fast, slow) -> -1/0/1

Return `1` if the `fast` and `slow` traces are positively correlated,
`-1` if they are negatively correlated, and `0` if their maximum
absolute cross-correlation is exactly 0.
"""
function _apparent_polarity(fast, slow)
    minval, maxval = extrema(DSP.xcorr(Seis.trace(fast), Seis.trace(slow)))
    max_xc = abs(minval) > abs(maxval) ? minval : maxval
    sign(max_xc)
end

end # module
