import RecipesBase

"""
    plot(s::SeisSplit.Result; max_samples=1000, antialias=true) -> ::Plots.Plot

Plot the results of shear wave splitting analysis from the `SeisSplit.splitting` function.

# Keyword arguments
- `antialias = true`: Control whether decimation of traces for plotting purposes
  uses antialising to avoid spurious signals.  `true` will avoid spurious aliasing
  signals, whilst `false` may be quicker.
- `max_samples = 1000`: Maximum number of samples to plot in the trace and particle
  motion subplots.  Windows or traces with more than this number are decimated to
  below this number to speed up plotting.
"""
plot

RecipesBase.@recipe function f(s::Result; max_samples=1000, antialias=true)
    grid --> false
    framestyle --> :box
    colorbar := false
    layout := 8
    title := hcat("Input", "Source pol. pre-corr", "Source pol. post-corr",
                  "Fast-slow pre-corr", "Fast-slow corr", "Small eigenvalue",
                  "PM pre-corr", "PM post-corr", "")
    titlefontsize --> 11
    legendfontsize --> 8
    fontfamily --> "Helvetica"
    size --> (600, 600)


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

    ## Input traces in original orientation with analysis window highlighted
    # Analysis window in colour behind traces
    # RecipesBase.@series begin
    #     subplot := 1
    #     w1, w2 =s.window_start, s.window_end
    #     []
    # end
    # Traces
    RecipesBase.@series begin
        subplot := 1
        label := [s.trace1.sta.cha s.trace2.sta.cha]
        trace1, trace2 = _decimate_to_max_samples.((s.trace1, s.trace2), max_samples,
            antialias)
        Seis.times(trace1), [Seis.trace(trace1) Seis.trace(trace2)]
    end
    # Window limits
    RecipesBase.@series begin
        subplot := 1
        seriestype := :vline
        label := ""
        linecolor := :red
        linewidth := 1
        [window...]
    end

    ## Rotated to spol
    RecipesBase.@series begin
        subplot := 2
        spol, spol_90 = _decimate_to_max_samples.(
            Seis.rotate_through(n_orig, e_orig, 180-s.spol),
            max_samples, antialias, Ref(window))
        label := ""
        linecolor := [:blue :red]
        Seis.times(spol), [Seis.trace(spol), Seis.trace(spol_90)]
    end

    ## Corrected and rotated to spol
    RecipesBase.@series begin
        subplot := 3
        label := ""
        linecolor := [:blue :red]
        n, e = _decimate_to_max_samples.((n_orig, e_orig), max_samples, antialias, Ref(window))
        SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
        fast, slow = Seis.rotate_through(n, e, 180-s.spol)
        Seis.times(fast), [Seis.trace(fast), Seis.trace(slow)]
    end

    ## Fast and slow before dt shift
    sfast, sslow = Seis.rotate_through!.(_decimate_to_max_samples.((n_orig, e_orig), max_samples,
        antialias, Ref(window))..., s.phi_best)
    # Flip slow wave if opposite polarity to make it easier to see if traces line up
    slow_polarity = let p = _apparent_polarity(sfast, sslow)
        # Account for zero polarity if traces not correlated at all
        p == 0 ? 1 : p
    end
    RecipesBase.@series begin
        subplot := 4
        label := ""
        legend := false
        linecolor := [:blue :red]
        xlims := window
        Seis.times(sfast), [Seis.trace(sfast), slow_polarity.*Seis.trace(sslow)]
    end

    ## Fast and slow after dt shift
    RecipesBase.@series begin
        subplot := 5
        label := ""
        legend := false
        linecolor := [:blue :red]
        xlims := window
        n, e = _decimate_to_max_samples.((n_orig, e_orig), max_samples, antialias, Ref(window))
        # Note we are taking the splitting off, hence -δt
        SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
        fast, slow = Seis.rotate_through(n, e, s.phi_best)
        [Seis.times(fast), Seis.times(slow)], [Seis.trace(fast), slow_polarity.*Seis.trace(slow)]
    end

    ## λ₂ surface
    # 10σ - 100σ contours
    RecipesBase.@series begin
        subplot := 6
        seriestype := :contour
        levels := 10:10:100
        linewidth := 1
        linecolor := :gray
        linealpha := 0.6
        s.dt, s.phi, s.lam2
    end

    # 1σ – 10σ contours
    RecipesBase.@series begin
        subplot := 6
        seriestype := :contour
        levels := 1:10
        linewidth := 1
        linecolor := :black
        s.dt, s.phi, s.lam2
    end
    # 1σ contour
    RecipesBase.@series begin
        subplot := 6
        seriestype := :contour
        levels := [1]
        linewidth := 3
        linecolor := :black
        s.dt, s.phi, s.lam2
    end
    # Best value
    RecipesBase.@series begin
        subplot := 6
        seriestype := :scatter
        xguide --> "dt / s"
        yguide --> "phi / deg"
        xlims := extrema(s.dt)
        ylims := extrema(s.phi)
        label := ""
        [s.dt_best], [s.phi_best]
    end

    ## Particle motion before correction
    amax = 1.2*maximum(x -> sqrt(x[1]^2 + x[2]^2), zip(Seis.trace(e_orig), Seis.trace(n_orig)))
    RecipesBase.@series begin
        subplot := 7
        aspect_ratio := :equal
        label := ""
        linecolor --> :black
        xguide := traces_n_e ? "East" : e_orig.sta.cha
        yguide := traces_n_e ? "North" : n_orig.sta.cha
        xticks --> nothing
        yticks --> nothing
        xlims := (-amax, amax)
        ylims := (-amax, amax)
        x, y = _decimate_to_max_samples.((e_orig, n_orig), max_samples, antialias,
            Ref(window))
        Seis.trace(x), Seis.trace(y)
    end

    ## Particle motion after correction
    RecipesBase.@series begin
        subplot := 8
        aspect_ratio := :equal
        label := ""
        linecolor --> :black
        xguide := traces_n_e ? "East" : e_orig.sta.cha
        yguide := traces_n_e ? "North" : n_orig.sta.cha
        xticks --> nothing
        yticks --> nothing
        xlims := (-amax, amax)
        ylims := (-amax, amax)
        n, e = deepcopy.((n_orig, e_orig))
        SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
        y, x = _decimate_to_max_samples.((n, e), max_samples, antialias, Ref(window))
        Seis.trace(x), Seis.trace(y)
    end

    ## Empty plot at bottom right
    # RecipesBase.@series begin
    #     subplot := 9
    #     framestyle := :empty
    #     nothing
    # end
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
