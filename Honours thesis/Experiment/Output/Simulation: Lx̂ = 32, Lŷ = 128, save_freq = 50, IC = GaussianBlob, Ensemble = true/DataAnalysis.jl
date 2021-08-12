#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

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
for i ∈ 2:data[1]["no_of_sims"]
    plot!(top_layer, t, second_moms[:, 1, i], label = "Ensemble member "*string(i))
end
plot!(top_layer, t, ens_second_mom[:, 1], label = "Ensemble average", line = (:dash, 3, :black))

## Bottom layer
bottom_layer = plot(t, second_moms[:, 2, 1], 
                label = "Ensemble member "*string(1), 
                xlabel = "t",
                ylabel = "σ²(t)",
                legend = :bottomright)
for i ∈ 2:data[1]["no_of_sims"]
    plot!(bottom_layer, t, second_moms[:, 2, i], label = "Ensemble member "*string(i))
end
plot!(bottom_layer, t, ens_second_mom[:, 2], label = "Ensemble average", line = (:dash, 3, :black))

##
plot(top_layer, bottom_layer, layout = (2, 1), size = (1200, 1200))