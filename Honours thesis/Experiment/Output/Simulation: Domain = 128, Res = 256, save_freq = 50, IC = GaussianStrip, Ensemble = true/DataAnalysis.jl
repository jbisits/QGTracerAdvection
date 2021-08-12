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
area_per = tracer_area_percentile(data; Cₚ = 0.5)
avg_area_per = avg_ensemble_tracer_area(data; Cₚ = 0.5)

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

diff = diffusivity(data, [21 70; 21 70]; Cₚ = 0.5)

ensemble_area_per = plot(upper_area, lower_area, layout = (2, 1), size = (1200, 1200))
save("tenensemble_area_per.png", ensemble_area_per)
##

## Ensemble area (Growth of average of concentration field)
t = time_vec(data[1])
area_per = tracer_area_percentile(data; Cₚ = 0.05)
ensemble_conc = ensemble_concentration(data)
ensemble_area = tracer_area_percentile(ensemble_conc; Cₚ = 0.05)

ensemble_upper = plot(t, ensemble_area[:, 1], 
                    title = "Growth of 50% of ensemble area of tracer patch in upper layer",
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
                    title = "Growth of 50% of ensemble area of tracer patch in lower layer",
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

ens_diff = diffusivity(data, [101 1; 101 1]; Cₚ = 0.05)

## Look at average area and second moment

second_moms = tracer_second_mom(data)
ens_conc = ensemble_concentration(data)
ens_second_mom = tracer_second_mom(ens_conc)

t = time_vec(data[1])

## Top layer
top_layer = plot(t, second_moms[:, 1, 1], 
                label = "Ensemble member "*string(1), 
                xlabel = "t",
                ylabel = "σ²(t)",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(top_layer, t, second_moms[:, 1, i], label = "Ensemble member "*string(i))
end
plot!(top_layer, t, ens_second_mom[:, 1], label = "Ensemble average", line = (:dash, 3, :black))

## Bottom layer
bottom_layer = plot(t, second_moms[:, 2, 1], 
                label = "Ensemble member "*string(1), 
                xlabel = "t",
                ylabel = "σ²(t)",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(bottom_layer, t, second_moms[:, 2, i], label = "Ensemble member "*string(i))
end
plot!(bottom_layer, t, ens_second_mom[:, 2], label = "Ensemble average", line = (:dash, 3, :black))

## 
plot(top_layer, bottom_layer, layout = (2, 1), size = (1200, 1200))