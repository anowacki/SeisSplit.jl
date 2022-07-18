# Quality control criteria for results

"""
    lam2_ratio(s::Result)

Return the ratio λ₂⁰/minimum(λ₂), where λ₂⁰ is the uncorrected λ₂.
This should be **maximised** when removing the best-fitting splitting operator,
turning elliptical particle motion linear.
"""
lam2_ratio(s::Result) = first(s.lam2)/minimum(s.lam2)

"""
    lam1_ratio(s::Result)

Return the ratio λ₁⁰/maximum(λ₁), where λ₁⁰ is the uncorrected λ₁.
This should be **minimised** when removing the best-fitting splitting operator,
turning elliptical particle motion linear.
"""
lam1_ratio(s::Result) = first(s.lam1)/maximum(s.lam1)

"""
    lam_ratio(s::Result)

Return the ratio (λ₂⁰/λ₁⁰) / [minimum(λ₂)/maximum(λ₁)].
This should be **maximised** when removing the best-fitting splitting operator,
turning elliptical partice motion linear.
"""
lam_ratio(s::Result) = lam2_ratio(s)/lam1_ratio(s)

"""
    quality(result) -> Q

Return the quality index `Q`, which is 1 for a perfect splitting measurement,
-1 for a perfect null, and 0 for a poor measurement.

See Wuestefeld et al. (Geophysical Prospecting, 2010) for details.

# Reference
- Wuestefeld, A., Al-Harrasi, O., Verdon, J.P., Wookey, J., Kendall, J.-M., 2010.
  A strategy for automated analysis of passive microseismic data to image
  seismic anisotropy and fracture characteristics.
  Geophysical Prospecting 58, 753–771.
  doi:[10.1111/j.1365-2478.2010.00891.x](https://doi.org/10.1111/j.1365-2478.2010.00891.x)
"""
function quality(s::Result)::Float64
    s.xcorr_dt === nothing &&
        throw(ArgumentError("rotation correlation map not computed for this result"))
    Δ = s.xcorr_dt_best/s.dt_best
    ϕdiff = abs(Seis.angle_difference(s.phi_best, s.xcorr_phi_best, true))
    ϕdiff > 45 && (ϕdiff = 90 - ϕdiff)
    Ω = ϕdiff/45
    d_null = sqrt(Δ^2 + (Ω - 1)^2)*√2
    d_null > 1 && (d_null = 1.0)
    d_good = sqrt((Δ - 1)^2 + Ω^2)*√2
    d_good > 1 && (d_good = 1.0)
    Q = d_null < d_good ? -(1 - d_null) : 1 - d_good
    dt_max = maximum(s.dt)
    Δ > 0.9 && s.dt_best/dt_max > 0.9 && (Q = -Q*(10 - 10*s.dt_best/dt_max))
    Q
end

"""
    snr_restivo_helffrich(result) -> SNR

Return the signal-to-noise ratio as defined by Restivo & Helffrich (GJI, 1999):
the signal is the maximum absolute amplitude of the corrected trace in the
source polarisation orientation, and the noise is the 2σ.

# Reference
- Restivo, A., Helffrich, G., 1999. Teleseismic shear wave splitting
  measurements in noisy environments. Geophys J Int 137, 821–830.
  doi:[10.1046/j.1365-246x.1999.00845.x](https://doi.org/10.1046/j.1365-246x.1999.00845.x)
"""
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
