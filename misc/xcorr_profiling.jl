# Try to recreate the example from Wuestefeld & Bokelmann, BSSA, 2007,
# Figure 1, and use this to both profile the splitting calculation
# overall, and also verify the rotation-correlation analysis.

using Revise, Seis, SeisSplit, Plots
import Juno, DelimitedFiles

# Synthetic parameters
noise = √(1/15)
spol = 60
freq = 1/8
phi_true = 0
dt_true = 1.3
dt_max = 3

# Create synthetics using programs from
#   https://github.com/anowacki/seismo-fortran
run(`create_wave -n $noise -f $freq -spol $spol`)
run(`splitwave $phi_true $dt_true`)

# Filter and save
t = read_sac("wave.BH[ENZ]", ".", echo=false)
bandpass!.(t, 0.02, 1, twopass=true)
write_sac.(t, t.meta.file)

# Calculate splitting using SHEBA for reference
open("sheba.in", "w") do f
    println(f, """
        # sheba.in
        wave
        BHE
        BHN
        BHZ
        1
        1 1
        $dt_max
        0
        0
        """)
end
withenv(()->run(`sheba_exec`), "SHEBA_WRITE_XC"=>"1")
sheba_phi, sheba_dphi, sheba_dt, sheba_ddt =
    parse.(Float64, split(readlines("wave_sheba.result")[2])[11:14])
sheba_Q = parse(Float64, split(readlines("wave_sheba.stats")[2])[5])
sheba_XC = DelimitedFiles.readdlm("sheba.xc")
sheba_phi_XC, sheba_dt_XC = parse.(Float64,
    split(readlines("sheba.xc_best")[2])[[1,3]])
sheba_phis_XC = -90:90
sheba_dts_XC = range(0, step=2t[1].delta, length=size(sheba_XC, 2))

# Calculate splitting using SeisSplit
e, n = Seis.read_sac("wave.BH[EN]", ".", echo=false)
@time s = splitting(n, e, e.picks.A.time, e.picks.F.time, dt_max=dt_max)

# Uncomment for profiling display with Juno
# Juno.@profiler begin
#     a = 0.0
#     for _ in 1:20
#         t, phi, xc, phi_best, dt_best = SeisSplit.rotation_correlation(n, e,
#             e.picks.A.time, e.picks.F.time)
#         a += t[1]
#     end
#     a
# end

p = plot(layout=(2,1))

# Xcorr heatmap, λ₂ contours, and best results
heatmap!(p[1], s.xcorr_dt, s.xcorr_phi, s.xcorr_map', c=:RdBu,
    xlabel="\\delta t / s", ylabel="\\phi / °",
    framestyle=:box, xlim=(0,dt_max), ylim=(-90,90))
contour!(p[1], s.dt, s.phi, s.lam2.*maximum(s.xcorr_map)./maximum(s.lam2))
scatter!(p[1], [s.xcorr_dt_best], [s.xcorr_phi_best], label="XC", markersize=15)
scatter!(p[1], [sheba_dt_XC], [sheba_phi_XC], label="SHEBA XC", markersize=12)
scatter!(p[1], [s.dt_best], [s.phi_best], label="EV", markersize=10)
scatter!(p[1], [sheba_dt], [sheba_phi], label="SHEBA EV", markersize=6)
scatter!(p[1], [dt_true], [phi_true], label="True", markersize=4)

# Xcorr surfaces, all normalised
XC_sheba_norm = sheba_XC ./ maximum(abs, sheba_XC)
XC_seissplit_norm = s.xcorr_map./maximum(abs, s.xcorr_map)
XC_splitlab_norm = DelimitedFiles.readdlm("splitlab.xc")
XC_splitlab_norm ./= maximum(abs, XC_splitlab_norm)
XC_splitlab_norm = reverse(XC_splitlab_norm, dims=1)
splitlab_phis = -90:89
splitlab_dts = range(0, dt_max, length=size(XC_splitlab_norm, 2))
iphi, idt = Tuple(argmin(XC_splitlab_norm))
splitlab_phi = splitlab_phis[iphi]
splitlab_dt = splitlab_dts[idt]

contour!(p[2], sheba_dts_XC, sheba_phis_XC, XC_sheba_norm,
    l=:red, levels=-1:0.2:1, label="SHEBA", xlim=(0,dt_max))
contour!(p[2], s.xcorr_dt, s.xcorr_phi, XC_seissplit_norm',
    l=:black, levels=-1:0.2:1, label="SeisSplit")
contour!(p[2], splitlab_dts, splitlab_phis, XC_splitlab_norm,
    l=:blue, levels=-1:0.2:1, label="Splitlab")
scatter!(p[2], [sheba_dt s.xcorr_dt_best splitlab_dt],
    [sheba_phi s.xcorr_phi_best splitlab_phi], markercolor=[:red :black :blue],
    markersize=[10 8 6],
    label=["SHEBA" "SeisSplit" "Splitlab"])

@show s.xcorr_phi_best, s.xcorr_dt_best
@show sheba_phi_XC, sheba_dt_XC
@show SeisSplit.quality(s), sheba_Q

p # show plot
