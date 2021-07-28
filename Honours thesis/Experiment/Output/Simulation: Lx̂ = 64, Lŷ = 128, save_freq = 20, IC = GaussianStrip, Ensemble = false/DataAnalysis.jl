#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 64, Lŷ = 128, save_freq = 20, IC = GaussianStrip, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

#Produce histogram plots from the saved concentration data
histplots = hist_plot(data; plot_freq = 1000, xlims_same = false)

#Produce heatmaps of tacer concentration from the saved concentration data
tracerplots = tracer_plot(data; plot_freq = 1000)

#Plot heatmaps and histograms togehter.
uppertacer = plot(tracerplots[1]...)

upperhist = plot(histplots[1]...)
plot(uppertacer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[2]...)

lowerhist = plot(histplots[2]...)
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

## Growth of 90% of tracer patch
t = time_vec(data)
area_per = tracer_area_percentile(data; conc_min = 0.1)
p1 = plot(t, area_per, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers; domain = 64 x 128. \n Gaussian strip IC",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of 90% of log(area of tracer patch) \n in both layers; domain = 64 x 128. \n Gaussian strip IC",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))