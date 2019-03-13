using Test
using Seis, SeisSplit
import Random, Dates

# Seed the random number generator 'randomly'
Random.seed!(Dates.second(Dates.now()))

"Test that all fields are approximately equal for two named tuples"
function splits_are_equal(s, s′)
    for k in keys(s)
        @test s[k] ≈ s′[k]
    end
end

@testset "All tests" begin
    @testset "Split" begin
        let t = Seis.read_sac.(joinpath(@__DIR__, "data", "wave.BH").*("E", "N"))
            e, n = t
            s = splitting(e, n)
            # Ensure trace ordering doesn't matter
            s′ = splitting(n, e)
            splits_are_equal(s, s′)
            @test s.phi_best ≈ 40 atol=step(s.phi)
            @test s.dt_best ≈ 1.4 atol=step(s.dt)
            @test s.spol ≈ 10 atol=1.0
            # Ensure starting orientation doesn't matter
            s′ = splitting(Seis.rotate_through(n, e, rand(1:360))...)
            splits_are_equal(s, s′)
        end
    end
    @testset "Errors" begin
        # TODO: Fix uncertainties and remove _broken when done.
        let t = Seis.read_sac.(joinpath(@__DIR__, "data", "wave.BH").*("E", "N"))
            e, n = t
            s = splitting(e, n)
            @test_broken s.dphi ≈ 0.75 atol=0.05
            @test_broken s.ddt ≈ 0.013 atol=0.01
            @test_broken s.dpol ≈ 0.35 atol=0.05
        end
    end
end