#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 64, Lŷ = 128, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

## Load in the data. This is an ensemble simulation so now have an array of dictionaries.
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
first_moms = first_moment(data)
first_mom_upper = plot(t, first_moms[:, 2, 1], label = "Member 1", legend = :bottomright)
for i ∈ 2:length(data)
    plot!(first_mom_upper, t, first_moms[:, 2, i], label = "Memeber "*string(i))
end
first_mom_upper

ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

plot!(first_mom_upper, t, ensemble_avg[:, 2], label = "Ensemble average", line = (:dash, 2, :black))

Δt = t[70] - t[30]
ΔA = ensemble_avg[70, 1] - ensemble_avg[30, 1]
K = ΔA / (4 * π * Δt)