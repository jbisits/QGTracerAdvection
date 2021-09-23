cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = QGPV, Ensemble = true"))
file = joinpath(pwd(), "SimulationData.jld2")

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

## First moment 
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

Δt = t[30] - t[1]
ΔA = ensemble_avg[30, :] .- ensemble_avg[1, :]
K = ΔA ./ (4 * π * Δt)

tp = tracer_plot(data[1])
plot(tp[:, 1]..., size = (1200, 1200))

## Second moment
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

ΔA² = ens_sec_mom[30, :] - ens_sec_mom[1, :]
Δt = t[30] - t[1]
Lₓ = data[1]["grid/Lx"]
K = ΔA² / (Lₓ^2 * 8 * Δt)