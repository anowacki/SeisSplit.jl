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
    layout := (2,3)
    
    # Input traces in N-E orientation
    n_orig, e_orig = Seis.rotate_through(s.trace1, s.trace2, -s.trace1.sta.azi)

    ## Input traces in original orientation
    title := ["Input" "Source pol. pre-corr" "Source pol. post-corr" "PM pre-corr" "PM post-corr" "\\lambda_2"]
    RecipesBase.@series begin
        subplot := 1
        label := [s.trace1.sta.cha s.trace2.sta.cha]
        Seis.times(s.trace1), [Seis.trace(s.trace1) Seis.trace(s.trace2)]
    end
    
    # RecipesBase.@series begin
    #     subplot := 2
    #     [s.trace2]
    # end
    #
    ## Rotated to spol
    RecipesBase.@series begin
        subplot := 2
        spol, spol_90 = Seis.rotate_through(n_orig, e_orig, 180-s.spol)
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
        SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
        fast, slow = Seis.rotate_through(n, e, 180-s.spol)
        Seis.times(fast), [Seis.trace(fast), Seis.trace(slow)]
    end

    ## Particle motion before correction
    amax = 1.4*maximum(x -> sqrt(x[1]^2 + x[1]^2), zip(Seis.trace(e_orig), Seis.trace(n_orig)))
    RecipesBase.@series begin
        subplot := 4
        aspect_ratio := :equal
        label := ""
        linecolor --> :black
        xlim := (-amax, amax)
        ylim := (-amax, amax)
        Seis.trace(e_orig), Seis.trace(n_orig)
    end

    ## Particle motion after correction
    RecipesBase.@series begin
        subplot := 5
        aspect_ratio := :equal
        label := ""
        linecolor --> :black
        xlim := (-amax, amax)
        ylim := (-amax, amax)
        n, e = deepcopy.((n_orig, e_orig))
        SeisSplit.apply_split!(n, e, s.phi_best, -s.dt_best)
        Seis.trace(e), Seis.trace(n)
    end

    ## λ₂ surface
    # 1σ – 10σ contours
    RecipesBase.@series begin
        subplot := 6
        seriestype := :contour
        levels --> 1:10
        linewidth := 1
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
        xlabel --> "\\delta t / s"
        ylabel --> "\\phi / °"
        xlim := extrema(s.dt)
        ylim := extrema(s.phi)
        seriestype := :scatter
        label := ""
        [s.dt_best], [s.phi_best]
    end
end
