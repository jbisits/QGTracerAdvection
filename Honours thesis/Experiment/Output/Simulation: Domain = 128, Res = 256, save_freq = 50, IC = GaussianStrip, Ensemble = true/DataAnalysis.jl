cd(joinpath(SimPath, "Output/Simulation: Domain = 128, res = 256, save_freq = 50, IC = GaussianStrip, Ensemble = true"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

t = time_vec(data)
area_per = tracer_area_percentile(data; conc_min = 0.1)
p1 = plot(t, area_per, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers; domain = 128, res = 256. \n Gaussian strip IC",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of 90% of log(area of tracer patch) \n in both layers; domain = 128, res = 256. \n Gaussian strip IC",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

file1 = joinpath(pwd(), "SimulationData_1.jld2")
data1 = load(file1)

t1 = time_vec(data1)
area_per1 = tracer_area_percentile(data1; conc_min = 0.1)
p11 = plot(t, area_per1, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers; domain = 128, res = 256. \n Gaussian strip IC",
        legend = :topleft
        )
logp11 = plot(t1, log.(area_per1), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of 90% of log(area of tracer patch) \n in both layers; domain = 128, res = 256. \n Gaussian strip IC",
            legend = :bottomright
            )
plot(p11, logp11, layout = (2, 1), size = (700, 700))
