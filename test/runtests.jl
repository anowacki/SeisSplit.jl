using Test
using Seis, SeisSplit
import Random, Dates

# Seed the random number generator 'randomly'
Random.seed!(Dates.second(Dates.now()))

read_test_data() = Seis.read_sac.(joinpath(@__DIR__, "data", "wave.BH").*("E", "N"))

"Test that all fields are approximately equal for two named tuples"
function splits_are_equal(s, s′)
    for f in fieldnames(SeisSplit.Result)
        f in (:trace1, :trace2) && continue # Don't compare traces themselves
        @test getfield(s, f) ≈ getfield(s′, f)
    end
end

@testset "All tests" begin
    @testset "Split" begin
        let t = read_test_data()
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
    @testset "Result" begin
        let t = read_test_data(), nphi=361, ndt=80, wb=10, we=35
            e, n = t
            s = splitting(e, n, 10, 35, nphi=nphi, ndt=ndt)
            @test length(s.phi) == nphi
            @test length(s.dt) == ndt
            @test size(s.lam1) == size(s.lam2) == (nphi, ndt)
            @test s.window_start == wb
            @test s.window_end == we
            @test s.trace1 == n
            @test s.trace2 == e
        end
    end
    @testset "Window" begin
        let t = read_test_data()
            e, n = t
            s = splitting(e, n, Seis.starttime(e), Seis.endtime(e))
            s′ = splitting(e, n)
            @test s.window_start == s′.window_start
            @test s.window_end == s′.window_end
        end
    end
    @testset "Errors" begin
        # TODO: Fix uncertainties and remove _broken when done.
        let t = read_test_data()
            e, n = t
            s = splitting(e, n, e.meta.SAC_user0, e.meta.SAC_user2)
            @test s.dphi ≈ 0.75 atol=0.5
            @test s.ddt ≈ 0.013 atol=0.1
            @test s.dspol ≈ 0.35 atol=0.5
        end
    end
end