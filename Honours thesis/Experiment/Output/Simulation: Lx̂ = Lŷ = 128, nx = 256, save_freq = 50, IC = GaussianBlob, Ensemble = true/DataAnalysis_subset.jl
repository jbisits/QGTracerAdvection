#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

## Load in the data for delay_time = Δt * 6000 with the new parameters, 
#seed = 1234 for first 10, seed = 4321 for 10-20, seed = 2341 for 20-30, seed = 3142 for 30-40, seed = 3241 for 40-50
data = Array{Dict{String, Any}}(undef, 50)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end
##
t = time_vec(data[1])
first_moms = first_moment(data)
ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

## Generate a subsets of the data over the grid and see how well area diagnostic performs

first_moments_subset = first_moment(ensemble_conc, 3 * 4, 4)

plot(t, first_moments_subset, 
    xlabel = "t",
    ylabel = "⟨A⟩",
    title = "Subset of ensemble data to compute ⟨A⟩",
    label = ["Upper layer" "Lower layer"],
    legend = :topleft)

plot(t, ensemble_avg, 
    xlabel = "t",
    ylabel = "⟨A⟩",
    title = "Ensemble data and subset of\nensemble data to compute ⟨A⟩",
    label = "Full data",
    legend = :topleft)
plot!(t, first_moments_subset, label = "Subset of data")

Δt = t[41] - t[1]
ΔA = first_moments_subset[41, :] .- first_moments_subset[1, :]
K_subset = ΔA ./ (4 * π * Δt)

K_sub_dim = @. K_subset * dims["Ld"] * 0.02

# Linear fits to the subset data

ens_fit_sub = [ones(length(t)) t] \ first_moments_subset

upperlinfit = plot(t, ens_fit_sub[1, 1] .+ ens_fit_sub[2, 1] .* t, 
                    title = "Upper layer linear best fit of the growth\n of subset ensemble average area",
                    label = "Best fit of subset ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(t, first_moments_subset[:, 1], 
    label = "Subset ensemble data")
#savefig(upperlinfit, "upperlinfitblob.png")

plot(t, ens_fit_sub[1, 2] .+ ens_fit_sub[2, 2] .* t, 
                    title = "Lower layer linear best fit of the growth\n of subset ensemble average area",
                    label = "Best fit of subset ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(t, first_moments_subset[:, 2], 
    label = "Subset ensemble data")
    
K_sub_dim = @. ( ens_fit_sub[2, :] / (4 * π) ) * dims["Ld"] * 0.02
