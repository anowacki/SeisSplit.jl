# Functions to which methods are added when loading a Makie
# backend and activating the SeisSplitMakieExt extension

"""
    plot_result(s::SeisSplit.Result; max_samples=1000, antialias=true) -> ::Makie.Figure
    Makie.plot(s::SeisSplit.Result; max_samples=1000, antialias=true) -> ::Makie.Figure

Plot the results of shear wave splitting analysis from the `SeisSplit.splitting` function.

# Keyword arguments
- `antialias = true`: Control whether decimation of traces for plotting purposes
  uses antialising to avoid spurious signals.  `true` will avoid spurious aliasing
  signals, whilst `false` may be quicker.
- `max_samples = 1000`: Maximum number of samples to plot in the trace and particle
  motion subplots.  Windows or traces with more than this number are decimated to
  below this number to speed up plotting.
"""
function plot_result end
