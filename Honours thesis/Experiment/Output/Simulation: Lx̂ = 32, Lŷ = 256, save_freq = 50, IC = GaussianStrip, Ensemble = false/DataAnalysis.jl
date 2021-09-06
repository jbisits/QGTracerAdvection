#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 256, save_freq = 50, IC = GaussianStrip, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

tracerplots = tracer_plot(data; plot_freq = 2000)
plot(tracerplots[:, 1]..., layout = (2, 3), size = (1200, 1200))

t = time_vec(data; days = true)
area_per = tracer_area_percentile(data; Cₚ = 0.5)

p1 = plot(t, area_per, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers; domain = 32 x 256. \n Gaussian strip IC",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of 90% of log(area of tracer patch) \n in both layers; domain = 32 x 256. \n Gaussian strip IC",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

## These are plots to compare over the different domain sizes and require code from the other sims, but the plots are saved
upperarea256 = plot(t, area_per[:, 1] .* 2^3,
                    title = "Growth of 90% of tracer area over four domains in \n upper layer where only meriodional length changes",
                    label = "32 x 256",
                    legend = :bottomright
                    )

lowerarea256 = plot(t, area_per[:, 2] .* 2^3,
                    title = "Growth of 90% of tracer area over four domains in \n lower layer where only meriodional length changes",
                    label = "32 x 256",
                    legend = :bottomright
                    )

save("upper_area4domains.png", upperarea256)
save("lower_area4domains.png", lowerarea256)

## Look at average area and second moment

t = time_vec(data)
first_mom = first_moment(data)
second_mom = second_moment(data)
plot(t, first_mom)
plot(t, second_mom)
## Look at concentration ~ grid cells

concVarea = concarea_animate(data)
mp4(concVarea, "ConcVsArea.mp4", fps = 18)