cd(joinpath(SimPath, "Output/Simulation: Domain = 128, res = 256, save_freq = 50, IC = GaussianStrip, Ensemble = true"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data. This is an ensemble simulation so now have an array of dictionaries.
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end

t = time_vec(data[1])
## Average of area percentage growth
area_per = tracer_area_percentile(data; conc_min = 0.1)
avg_area_per = avg_ensemble_tracer_area(data; conc_min = 0.1)

upper_area = plot(t, area_per[:, 1, 1], 
                label = "Ensemble member 1", 
                title = "Growth of 90% of area of tracer patch (upper layer) \n from ensemble simulation with 10 memebers",
                legend = :topleft)
for i ∈ 2:length(data)
    plot!(upper_area, t, area_per[:, 1, i], 
        label = "Ensemble member "*string(i)
        )
end
plot!(upper_area, t, avg_area_per[:, 1], line = (:dash, 2, :black), label = "Average")

lower_area = plot(t, area_per[:, 2, 1], 
                label = "Ensemble member 1", 
                title = "Growth of 90% of area of tracer patch (lower layer) \n from ensemble simulation with 10 memebers",
                legend = :topleft)
for i ∈ 2:length(data)
    plot!(lower_area, t, area_per[:, 2, i], 
        label = "Ensemble member "*string(i)
        )
end  
plot!(lower_area, t, avg_area_per[:, 2], line = (:dash, 2, :black), label = "Average")

diff = diffusivity(data, [21 70; 21 70]; conc_min = 0.1)

ensemble_area_per = plot(upper_area, lower_area, layout = (2, 1), size = (1200, 1200))
save("tenensemble_area_per.png", ensemble_area_per)
##

## Ensemble area (Growth of average of concentration field)
t = time_vec(data[1])
area_per = tracer_area_percentile(data; conc_min = 0.1)
ensemble_area = ensemble_average_area(data; conc_min = 0.1)

ensemble_upper = plot(t, ensemble_area[:, 1], 
                    title = "Growth of 90% of ensemble area of tracer patch in upper layer",
                    legend = :topleft,
                    label = "Ensemble area",
                    line = (:dash, 2, :black)
                    )
for i ∈ 1:length(data)
    plot!(ensemble_upper, t, area_per[:, 1, i], 
        label = "Ensemble member "*string(i)
        )
end
ensemble_lower = plot(t, ensemble_area[:, 2], 
                    title = "Growth of 90% of ensemble area of tracer patch in lower layer",
                    legend = :topleft,
                    label = "Ensemble area",
                    line = (:dash, 2, :black)
                    )
for i ∈ 1:length(data)
    plot!(ensemble_lower, t, area_per[:, 2, i], 
        label = "Ensemble member "*string(i)
        )
end
ens_plot = plot(ensemble_upper, ensemble_lower, layout = (2,1), size = (1200, 1200))
save("ensembleavg.png", ens_plot)