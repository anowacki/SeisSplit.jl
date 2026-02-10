@testset "QC" begin
    let t = read_test_data()
        e, n = t
        s = splitting(e, n, e.meta.SAC_user0, e.meta.SAC_user2)
        s′ = splitting(Seis.rotate_through(e, n, rand(1:360))...,
                       e.meta.SAC_user0, e.meta.SAC_user2)

        @testset "Restivo & Helffrich SNR" begin
            @test snr_restivo_helffrich(s) ≈ 19.695 atol=2
            @test snr_restivo_helffrich(s) ≈ snr_restivo_helffrich(s′)
        end

        @testset "Q" begin
            @test quality(s) ≈ 0.902 atol=0.02
            @test quality(s) ≈ quality(s′)
        end

        @testset "Q error without xcorr" begin
            s″ = splitting(n, e, xcorr=false)
            @test_throws ArgumentError quality(s″)
        end
    end
end
