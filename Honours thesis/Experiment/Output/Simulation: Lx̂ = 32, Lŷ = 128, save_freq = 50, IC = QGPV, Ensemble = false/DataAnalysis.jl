cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = QGPV, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

## Look at first moment.
t = time_vec(data)
first_mom = first_moment(data)
plot(t, first_mom, label = false)

## And the second
second_mom = second_moment(data)
plot(t, second_mom, label = false)

## First QGPV experiment so look at what it looks like
QGPV_plot = tracer_plot(data)
upperplots = plot(QGPV_plot[:, 1]...)
lowerplots = plot(QGPV_plot[:, 2]...)
plot(upperplots, lowerplots, layout = (2, 1), size = (1200, 1200))