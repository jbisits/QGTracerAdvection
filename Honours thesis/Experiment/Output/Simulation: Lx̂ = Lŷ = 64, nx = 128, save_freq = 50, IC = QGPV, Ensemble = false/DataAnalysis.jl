cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 64, nx = 128, save_freq = 50, IC = QGPV, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

PVplots = tracer_plot(data; plot_freq = 1000)
plot(PVplots[:, 1]..., size = (1200, 1200))
plot(PVplots[:, 2]..., size = (1200, 1200))