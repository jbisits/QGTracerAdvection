#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = GaussianStrip, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

t = time_vec(data)
mer_sec_mom = meridional_second_mom(data)

plot(t, mer_sec_mom, label = false)


## Older probably not sueful stuff
area_per = tracer_area_percentile(data; Cₚ = 0.1)
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

plot!(upperarea256, t, area_per[:, 1] .* 2^2, label = "32 x 128")
plot!(lowerarea256, t, area_per[:, 2] .* 2^2, label = "32 x 128")

## Look at average area and second moment

tsecs = time_vec(data; time_measure = "secs")
tdays = time_vec(data; time_measure = "days")
avg_area = tracer_avg_area(data)
second_moments = tracer_second_mom(data)

plot(tdays, avg_area, label = false)
plot(tdays, second_moments ,label = false)

K = second_moments[2:end, :] ./ (2 .* tsecs[2:end])
plot(tdays[2:end], K, label = false)