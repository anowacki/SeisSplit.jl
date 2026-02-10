# Script to calculate the relative speed of computing the cross correlation
# between two vectors in the time domain (using SeisSplit.cross_correlation)
# and frequency domain (using DSP.xcorr).
#
# This is useful so that we can call the appropriate type of cross-
# correlation calculation for the rotation-correlation approach in order to
# avoid the FFT-iFFT overhead in always using a frequency-domain approach.

# Results on 2019-03-025 on a MacBook Pro 2.9 GHz Intel Core i5:
#
# Float32
# 300: DSP   135.636 μs (144 allocations: 42.20 KiB)
# 300: SeisSplit:   78.421 μs (0 allocations: 0 bytes)
# 350: DSP   136.066 μs (144 allocations: 42.33 KiB)
# 350: SeisSplit:   107.384 μs (0 allocations: 0 bytes)
# 400: DSP   136.238 μs (144 allocations: 42.52 KiB)
# 400: SeisSplit:   141.003 μs (0 allocations: 0 bytes)
# 450: DSP   136.002 μs (144 allocations: 42.80 KiB)
# 450: SeisSplit:   179.090 μs (0 allocations: 0 bytes)
# 500: DSP   136.445 μs (144 allocations: 42.94 KiB)
# 500: SeisSplit:   221.634 μs (0 allocations: 0 bytes)
#
# Float64
# 300: DSP   161.177 μs (145 allocations: 75.34 KiB)
# 300: SeisSplit:   78.634 μs (0 allocations: 0 bytes)
# 350: DSP   161.443 μs (145 allocations: 75.78 KiB)
# 350: SeisSplit:   107.701 μs (0 allocations: 0 bytes)
# 400: DSP   160.819 μs (145 allocations: 76.03 KiB)
# 400: SeisSplit:   141.280 μs (0 allocations: 0 bytes)
# 450: DSP   161.884 μs (145 allocations: 76.53 KiB)
# 450: SeisSplit:   179.604 μs (0 allocations: 0 bytes)
# 500: DSP   161.987 μs (145 allocations: 76.91 KiB)
# 500: SeisSplit:   222.271 μs (0 allocations: 0 bytes)


import SeisSplit, DSP
using BenchmarkTools

times = Dict()

for T in (Float32, Float64)
    println(T)
    for n in 300:50:500
        print(n, ": DSP:       ")
        @btime c = DSP.xcorr(a, b) setup=(a=rand($T, $n); b=rand($T, $n); c=zeros($T, 2*$n-1)) seconds=1
        print(n, ": SeisSplit: ")
        @btime SeisSplit.cross_correlation!(c, a, b) setup=(a=rand($T, $n); b=rand($T, $n); c=zeros($T, 2*$n-1)) seconds=1
        a = rand(T, n)
        b = rand(T, n)
        @assert DSP.xcorr(a, b) ≈ SeisSplit.cross_correlation(a, b)
    end
    println()
end