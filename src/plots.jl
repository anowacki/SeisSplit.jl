import RecipesBase

"""
    plot(s::SeisSplit.Result) -> ::Plots.Plot

Plot the results of shear wave splitting analysis from the `SeisSplit.splitting` function.
"""
plot

RecipesBase.@recipe function f(s::Result)
    grid --> false
    framestyle --> :box
    colorbar := false
    layout := 8
    title := hcat("Input", "Source pol. pre-corr", "Source pol. post-corr",
                  "Fast-slow pre-corr", "Fast-slow corr", "\$\\lambda_2\$",
                  "PM pre-corr", "PM post-corr", "")
    titlefontsize --> 11
    legendfontsize --> 8

    # Input traces in N-E orientation
    n_orig, e_orig = Seis.rotate_through(s.trace1, s.trace2, -s.trace1.sta.azi)

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
        Seis.times(s.trace1), [Seis.trace(s.trace1) Seis.trace(s.trace2)]
    end
    # Window limits
    RecipesBase.@series begin
        subplot := 1
        seriestype := :vline
        label := ""
        linecolor := :red
        linewidth := 1
        [s.window_start, s.window_end]
    end

    ## Rotated to spol
    RecipesBase.@series begin
        subplot := 2
        spol, spol_90 = Seis.cut.(Seis.rotate_through(n_orig, e_orig, 180-s.spol),
                                 s.window_start, s.window_end)
        label := ""
        linecolor := [:blue :red]
        Seis.times(spol), [Seis.trace(spol), Seis.trace(spol_90)]
    end

    ## Corrected and rotated to spol
    RecipesBase.@series begin
        subplot := 3
        label := ""
        linecolor := [:blue :red]
        n, e = deepcopy.((n_orig, e_orig))
        SeisSplit.apply_split!(Seis.cut!.((n, e), s.window_start, s.window_end)...,
                               s.phi_best, -s.dt_best)
        fast, slow = Seis.rotate_through(n, e, 180-s.spol)
        Seis.times(fast), [Seis.trace(fast), Seis.trace(slow)]
    end

    ## Fast and slow before dt shift
    sfast, sslow = Seis.cut.(Seis.rotate_through(n_orig, e_orig, s.phi_best),
                             s.window_start, s.window_end)
    RecipesBase.@series begin
        subplot := 4
        label := ""
        legend := false
        linecolor := [:blue :red]
        Seis.times(sfast), [Seis.trace(sfast), Seis.trace(sslow)]
    end

    ## Fast and slow after dt shift
    RecipesBase.@series begin
        subplot := 5
        legend := false
        linecolor := [:blue :red]
        xlims := (s.window_start, s.window_end)
        # Fake the shift by changing sample times
        [Seis.times(sfast), Seis.times(sslow) .- s.dt_best], [Seis.trace(sfast), Seis.trace(sslow)]
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
        xguide --> "\$\\delta \\mathrm{t} \\, / \\, \\mathrm{s}\$"
        yguide --> "\$\\phi \\, / \\, ^\\circ\$"
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
        xguide := "East"
        yguide := "North"
        xticks --> nothing
        yticks --> nothing
        xlims := (-amax, amax)
        ylims := (-amax, amax)
        Seis.trace(Seis.cut(e_orig, s.window_start, s.window_end)),
            Seis.trace(Seis.cut(n_orig, s.window_start, s.window_end))
    end

    ## Particle motion after correction
    RecipesBase.@series begin
        subplot := 8
        aspect_ratio := :equal
        label := ""
        linecolor --> :black
        xguide := "East"
        yguide := "North"
        xticks --> nothing
        yticks --> nothing
        xlims := (-amax, amax)
        ylims := (-amax, amax)
        n, e = deepcopy.((n_orig, e_orig))
        SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
        Seis.cut!.((n, e), s.window_start, s.window_end)
        Seis.trace(e), Seis.trace(n)
    end

    ## Empty plot at bottom right
    # RecipesBase.@series begin
    #     subplot := 9
    #     framestyle := :empty
    #     nothing
    # end
end
