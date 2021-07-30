# This expermient had a flow with Ld = 2 in the simualtion, just interested how the mixing goes when the flow has larger eddies

#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 25, IC = GaussianStrip, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

t = time_vec(data)
area_per = tracer_area_percentile(data; conc_min = 0.1)
p1 = plot(t, area_per, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers; domain = 32 x 128. \n Gaussian strip IC",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of 90% of log(area of tracer patch) \n in both layers; domain = 32 x 128. \n Gaussian strip IC",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

diff = diffusivity(data, [85 285;1 2]; conc_min = 0.1)