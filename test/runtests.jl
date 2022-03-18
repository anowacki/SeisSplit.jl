using Test
using Seis, SeisSplit
import Random, Dates

# Seed the random number generator 'randomly'
Random.seed!(Dates.second(Dates.now()))

read_test_data() = Seis.read_sac.(joinpath(@__DIR__, "data", "wave.BH").*("E", "N"))

"Test that all fields are approximately equal for two `SeisSplit.Result`s"
function splits_are_equal(s, s′, atol=0.01)
    for f in fieldnames(SeisSplit.Result)
        f in (:trace1, :trace2, :xcorr_map, :reference_frame) && continue # Don't compare traces themselves
        @test getfield(s, f) ≈ getfield(s′, f) atol=atol
    end
    @test s.reference_frame == s′.reference_frame
end


@testset "All tests" begin
    include("minimum_eigenvalue.jl")
    include("rotation_correlation.jl")
    include("qc.jl")
end
