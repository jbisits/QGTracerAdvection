#This script is to check that the initial conditions are being set correctly and behaving when being advected.

## Gaussian strip IC
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 64, save_freq = 5, IC = GaussianStrip, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

#Produce histogram plots from the saved concentration data
histplots = hist_plot(data; plot_freq = 10, xlims_same = false)

#Produce heatmaps of tacer concentration from the saved concentration data
tracerplots = tracer_plot(data; plot_freq = 10)

#Plot heatmaps and histograms togehter.
uppertacer = plot(tracerplots[1]...)

upperhist = plot(histplots[1]...)
plot(uppertacer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[2]...)

lowerhist = plot(histplots[2]...)
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

## Gaussain blob IC
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 64, save_freq = 5, IC = GaussianBlob, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)
tracerplots = tracer_plot(data; plot_freq = 10)