import DelimitedFiles, DSP

@testset "Rotation correlation" begin
    @testset "Correlation" begin
        let Ts = (Float32, Float64), N = 1000
            for T in Ts
                a, b = rand(T, N), rand(T, N)
                xc = similar(a, 2N-1)
                work_arrays = SeisSplit._correlation_arrays(T, N)
                @test DSP.xcorr(a, b) ≈ SeisSplit._xcorr_fd!(xc, a, b, work_arrays...)
            end
        end
    end

    @testset "Xcorr" begin
        let t = read_test_data(), nphi = 180
            e, n = t
            win1, win2 = n.meta.SAC_user0, n.meta.SAC_user2
            dts, phis, xcorr_map, phi_best, dt_best =
                SeisSplit.rotation_correlation(n, e, win1, win2, nphi=nphi÷2)
            s = splitting(n, e, win1, win2, nphi=nphi)
            @test s.xcorr_phi == phis
            @test s.xcorr_dt == dts
            @test s.xcorr_map == xcorr_map
            @test s.xcorr_phi_best == phi_best
            @test s.xcorr_dt_best == dt_best
            @test phi_best ≈ 42 atol=1
            @test dt_best ≈ 1.35 atol=step(s.xcorr_phi)
        end
    end

    @testset "Xcorr map" begin
        # Compare to SHEBA-computed answer, though SHEBA appears
        # to include an extra point before 0 s dt, and takes the
        # absolute correlation coefficient.
        let t = read_test_data(), nphi = 180
            e, n = t
            dts, phis, xc, phi_best, dt_best =
                rotation_correlation(n, e, n.meta.SAC_user0, n.meta.SAC_user2, nphi=nphi)
            sheba_map = DelimitedFiles.readdlm(joinpath(@__DIR__, "data", "sheba.xc"))
            sheba_map = sheba_map'[2:end,1:180]
            xc = xc[1:2:80,:] .|> abs
            @test sheba_map ≈ xc rtol=0.01
        end
    end
end
