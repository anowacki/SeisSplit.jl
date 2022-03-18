"Maximum string length of the field names of `Result`"
const MAX_FIELD_NAME_LENGTH = maximum(x->length(String(x)), fieldnames(Result))

function Base.show(io::IO, ::MIME"text/plain", s::Result{T,V}) where {T,V}
    println("SeisSplit.Result{$T,$V}:")
    maxlen = MAX_FIELD_NAME_LENGTH
    for f in fieldnames(Result)
        unit = f in (:phi_best, :dphi, :spol, :dspol, :xcorr_phi_best) ? " °" :
               f in (:dt_best, :ddt, :xcorr_dt_best, :window_start, :window_end) ? " s" : ""
        if f in (:trace1, :trace2)
            t = getfield(s, f)
            out = "Seis.Trace($(Seis.channel_code(t)): " *
                "delta=$(t.delta), b=$(t.b), nsamples=$(Seis.nsamples(t)))"
            println(lpad(String(f), maxlen), ": ", out)
        elseif f in (:lam1, :lam2, :xcorr_dt, :xcorr_phi, :xcorr_map)
            λ = getfield(s, f)
            if !isnothing(λ)
                siz = join(size(λ), "×")
                println(io, lpad(String(f), maxlen), ": ", siz,
                    " Array{$(eltype(λ)),$(ndims(λ))}: [", λ[1], ", ..., ", λ[2], "]")
            else
                println(io, lpad(String(f), maxlen), ": (not calculated)")
            end
        else
            val = getfield(s, f)
            if !isnothing(val)
                println(io, lpad(String(f), maxlen), ": ", getfield(s, f), unit)
            else
                println(io, lpad(String(f), maxlen), ": (not calculated)")
            end
        end
    end
end

"""
    write_fasst_result_file(file, s::Result; header=trye)

Write the a shear wave splitting result `s` to `file` in FASST format.
If `header` is `false`, then do not write the header line.
"""
write_fasst_result_file(file, s::Result; header=true) =
    write(file, fasst_result_string(s; header=header))

"""
    fasst_result_string(s::Result; header=false) -> str

Return a string `str` in the format used the by the FASST shear wave splitting program
for a result `s`.  If `header` is `false`, then do not include the header line.
"""
function fasst_result_string(s::Result; header=true)
    _missing(s) = ismissing(s) ? "?" : s
    t = s.trace1
    fasst_string = ""
    if header
        fasst_string *= "#  NETWK  STATION  EVENT_LON  " *
                        "EVENT_LAT  EVT_DEP   STA_LON    STA_LAT DISTANCE AZIMUTH BACKAZI    " *
                        "FAST   DFAST    TLAG   DTLAG    SPOL   DSPOL  BEGIN_WIN    " *
                        "END_WIN   PHI_1   PHI_2 NPHI MAX_DT  NDT FILE1           FILE2\n"
    end
    # @printf cannot cope with splitting the format string across lines
    fasst_string *= Printf.@sprintf(
        "%8s %8s %10.5f %10.5f %7.2f %10.5f %10.5f %8.2f %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f %7.2f %7.2f %10.3f %10.3f %7.2f %7.2f %4i %6.2f %4i \"%s\" \"%s\"\n",
        _missing(t.sta.net), _missing(t.sta.sta), t.evt.lon, t.evt.lat,
        t.evt.dep, t.sta.lon, t.sta.lat, Seis.distance_deg(t),
        Seis.azimuth(t), Seis.backazimuth(t), s.phi_best, s.dphi,
        s.dt_best, s.ddt, s.spol, s.dspol, s.window_start,
        s.window_end, first(s.phi), last(s.phi), length(s.phi),
        maximum(s.dt), length(s.dt), s.trace1.meta.file, s.trace2.meta.file)
    fasst_string
end

"""
    write_fasst_lam2_file(file, s::Result)
    
Write the λ₂ grid for the splitting result `s` to `file` as a FASST-formatted file.
"""
function write_fasst_lam2_file(file, s::Result)
    open(file, "w") do f
        DelimitedFiles.writedlm(f, hcat(length(s.phi), length(s.dt),
                                        first(s.phi), last(s.phi),
                                        last(s.dt), minimum(s.lam2), maximum(s.lam2)), ' ')
        DelimitedFiles.writedlm(f, step(s.dt))
        DelimitedFiles.writedlm(f, hcat(-1.0, Seis.azimuth(s.trace1),
            Seis.backazimuth(s.trace1)), ' ')
        DelimitedFiles.writedlm(f, s.lam2, ' ')
    end
end
