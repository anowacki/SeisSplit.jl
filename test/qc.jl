@testset "QC" begin
    let t = read_test_data()
        e, n = t
        s = splitting(e, n, e.meta.SAC_user0, e.meta.SAC_user2)
        @test snr_restivo_helffrich(s) ≈ 19.695 atol=2
        s′ = splitting(Seis.rotate_through(n, e, rand(1:360))...,
                       e.meta.SAC_user0, e.meta.SAC_user2)
        @test snr_restivo_helffrich(s) ≈ snr_restivo_helffrich(s′)
        @test quality(s) ≈ 0.902 atol=0.02
        @test quality(s) ≈ quality(s′)
        s″ = splitting(n, e, xcorr=false)
        @test_throws ArgumentError quality(s″)
    end
end
