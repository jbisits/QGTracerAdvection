cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = GaussianStrip, Ensemble = true"))

## Load in the data. This is an ensemble simulation so now have an array of dictionaries.

#The saved data `Simulationdata` to `Simulationdata_9` are one set of ensemble expt with delay time = Δt * 3000
#The saved data `Simulationdata_10` to `Simulationdata_19` are one set of ensemble expt with delay time = Δt * 5000 to see if more developed flow changes anything.
data = Array{Dict{String, Any}}(undef, 10)
#Read in first ensemble sim
for i ∈ 1:length(data)
    if i == 1 
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end
#Read in second ensemble sim
for i ∈ 1:length(data)
    if i == 1 
        file = joinpath(pwd(), "SimulationData_10.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 9)*".jld2")
        data[i] = load(file)
    end
end

##
t = time_vec(data[1])
sec_mom = second_moment(data)

upperplot = plot(t, sec_mom[:, 1, 1],
                title = "Upper layer",
                xlabel = "t",
                ylabel = "σ²ₐ",
                label = "Member 1",
                legend = :bottomright)
lowerplot = plot(t, sec_mom[:, 2, 1],
                title = "Lower layer",
                xlabel = "t",
                ylabel = "σ²ₐ",
                label = "Member 1",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(upperplot, t, sec_mom[:, 1, i], label = "Member "*string(i)) 
    plot!(lowerplot, t, sec_mom[:, 2, i], label = "Member "*string(i)) 
end

ens_conc = ensemble_concentration(data)
ens_sec_mom = second_moment(ens_conc)

plot!(upperplot, t, ens_sec_mom[:, 1], label = "Ensemble", line = (:dash, :black, 2))
plot!(lowerplot, t, ens_sec_mom[:, 2], label = "Ensemble", line = (:dash, :black, 2))

plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))

ΔA² = ens_sec_mom[50, :] - ens_sec_mom[10, :]
Δt = t[50] - t[10]
Lₓ = data[10]["grid/Lx"]
K = ΔA² / (Lₓ^2 * 8 * Δt)
##
tp = tracer_plot(data[1])
plot(tp[:, 1]..., size = (1200, 1200))
##
t = time_vec(data[1])
mer_sec_mom = meridional_second_mom(data)

upperplot = plot(t, mer_sec_mom[:, 1, 1],
                title = "Upper layer",
                xlabel = "t",
                ylabel = "σ²y",
                label = "Member 1",
                legend = :bottomright)
lowerplot = plot(t, mer_sec_mom[:, 2, 1],
                title = "Lower layer",
                xlabel = "t",
                ylabel = "σ²y",
                label = "Member 1",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(upperplot, t, mer_sec_mom[:, 1, i], label = "Member "*string(i)) 
    plot!(lowerplot, t, mer_sec_mom[:, 2, i], label = "Member "*string(i)) 
end

ens_conc = ensemble_concentration(data)
ens_mer_sec_mom = meridional_second_mom(ens_conc)

plot!(upperplot, t, ens_mer_sec_mom[:, 1], label = "Ensemble", line = (:dash, :black, 2))
plot!(lowerplot, t, ens_mer_sec_mom[:, 2], label = "Ensemble", line = (:dash, :black, 2))

plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))
