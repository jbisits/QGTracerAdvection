#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

#The saved data `Simulationdata` to `Simulationdata_9` are ensemble with delay_time = Δt * 3000 
#The saved data `Simulationdata_10` to `Simulationdata_19` are ensemble with delay_time = Δt * 5000 
#The saved data `Simulationdata_20` to `Simulationdata_29` are ensemble with delay_time = Δt * 4000 
##Read in the first ensemble sim
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
## Read in second ensemble sim
for i ∈ 1:length(data)
    if i == 1 
        file = joinpath(pwd(), "SimulationData_10.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 9)*".jld2")
        data[i] = load(file)
    end
end
## Read in third ensemble sim
for i ∈ 1:length(data)
    if i == 1 
        file = joinpath(pwd(), "SimulationData_20.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 19)*".jld2")
        data[i] = load(file)
    end
end

##
t = time_vec(data[1])
first_moms = first_moment(data)
first_mom_upper = plot(t, first_moms[:, 1, 1], 
                        xlabel = "t",
                        ylabel  = "⟨A⟩",
                        title = "Average area growth in upper layer",
                        label = "Member 1",
                        legend = :bottomright)
first_mom_lower = plot(t, first_moms[:, 2, 1], 
                        xlabel = "t",
                        ylabel  = "⟨A⟩",
                        title = "Average area growth in lower layer",
                        label = "Member 1", 
                        legend = :bottomright)
for i ∈ 2:length(data)
    plot!(first_mom_upper, t, first_moms[:, 1, i], label = "Memeber "*string(i))
    plot!(first_mom_lower, t, first_moms[:, 2, i], label = "Memeber "*string(i))
end

ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

plot!(first_mom_upper, t, ensemble_avg[:, 1], label = "Ensemble average", line = (:dash, 2, :black))
plot!(first_mom_lower, t, ensemble_avg[:, 2], label = "Ensemble average", line = (:dash, 2, :black))

plot(first_mom_upper, first_mom_lower, layout = (2, 1), size = (800, 800))
##
Δt = t[80] - t[1]
ΔA = ensemble_avg[80, :] .- ensemble_avg[1, :]
K = ΔA ./ (4 * π * Δt)

#Now depends on the value of U the background horizontal velocity but looking quite reasonable

##
tp = tracer_plot(data[1]; plot_freq = 1000)
plot(tp[:, 1]..., size = (1200, 1200))