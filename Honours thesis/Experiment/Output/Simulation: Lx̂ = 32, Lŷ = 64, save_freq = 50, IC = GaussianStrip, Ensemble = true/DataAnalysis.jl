cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 64, save_freq = 50, IC = GaussianStrip, Ensemble = true"))

## Load in the data delat_time = Δt * 3000
data = Array{Dict{String, Any}}(undef, 20)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end

## Load in the data delay_time = Δt * 3500
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData_20.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 19)*".jld2")
        data[i] = load(file)
    end
end
## Load in the data delay_time = Δt * 4000
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData_30.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 29)*".jld2")
        data[i] = load(file)
    end
end
## Load in the data delay_time = Δt * 4500
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData_40.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 39)*".jld2")
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
                label = "Ensemble member",
                legend = :bottomright)
lowerplot = plot(t, sec_mom[:, 2, 1],
                title = "Lower layer",
                xlabel = "t",
                ylabel = "σ²ₐ",
                label = "Ensemble member",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(upperplot, t, sec_mom[:, 1, i], label = false) 
    plot!(lowerplot, t, sec_mom[:, 2, i], label = false) 
end

ens_conc = ensemble_concentration(data)
ens_sec_mom = second_moment(ens_conc)

plot!(upperplot, t, ens_sec_mom[:, 1], label = "Ensemble", line = (:dash, :black, 2))
plot!(lowerplot, t, ens_sec_mom[:, 2], label = "Ensemble", line = (:dash, :black, 2))

plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))

##
ΔA² = ens_sec_mom[25, :] - ens_sec_mom[1, :]
Δt = t[25] - t[1]
Lₓ = data[10]["grid/Lx"]
K = ΔA² / (Lₓ^2 * 8 * Δt)