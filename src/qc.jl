# Quality control criteria for results

"Return the ratio λ₂⁰/minimum(λ₂), where λ₂⁰ is the uncorrected λ₂.
This should be maximised when removing the best-fitting splitting operator
turns particle motion from elliptical to linear."
lam2_ratio(s::Result) = first(s.lam2)/minimum(s.lam2)

"Return the ratio λ₁⁰/maximum(λ₁), where λ₁⁰ is the uncorrected λ₁.
This should be minimised when removing the best-fitting splitting operator
turns particle motion from elliptical to linear."
lam1_ratio(s::Result) = first(s.lam1)/maximum(s.lam1)

"Return the ratio (λ₂⁰/λ₁⁰) / [minimum(λ₂)/maximum(λ₁)].
This should be maximised when removing the best-fitting splitting operator
turns partice motion from elliptical to linear."
lam_ratio(s::Result) = lam2_ratio(s)/lam1_ratio(s)

"Return the signal-to-noise ratio as defined by Restivo & Helffrich (GJI, 1999):
the signal is the maximum absolute amplitude of the corrected trace in the
source polarisation orientation, and the noise is the 2σ "
function snr_restivo_helffrich(s::Result)
    N, E = deepcopy.((s.trace1, s.trace2))
    Seis.cut!.((N, E), s.window_start, s.window_end)
    rotate_traces!(N, E, -N.sta.azi) # Now in N-E frame
    apply_split!(N, E, s.phi_best, -s.dt_best)
    rotate_traces!(N, E, -s.spol) # Now in spol-(spol+90) frame
    signal = maximum(abs, Seis.trace(N))
    noise = 2*Statistics.std(Seis.trace(E))
    signal/noise
end