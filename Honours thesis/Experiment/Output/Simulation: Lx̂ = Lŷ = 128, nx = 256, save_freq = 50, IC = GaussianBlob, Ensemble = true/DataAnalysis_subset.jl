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
ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

ϕ = π / 3
long_degrees = cos(ϕ) * 111.32
one_degree = 4 # 4 grid cells ≈ 60km which is an overestimate of 1 degree at this latitude (≈ 55km) but what I am using
## Generate a subsets of the data over the grid and see how well area diagnostic performs
no_of_degrees = 4

ens_avg_area_subset = first_moment(ensemble_conc, no_of_degrees * one_degree, 4)

plot(t, ens_avg_area_subset, 
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
plot!(t, ens_avg_area_subset, label = "Subset of data")

Δt = t[end] - t[round(Int64, 3*end / 4)]
ΔA = ens_avg_area_subset[end, :] .- ens_avg_area_subset[round(Int64, 3*end / 4), :]
K_subset = ΔA ./ (4 * π * Δt)

dims = nondim2dim(data[1])
K_sub_dim = @. K_subset * dims["Ld"] * 0.02

## Linear fits to the subset data

ens_fit_sub = [ones(length(t)) t] \ ens_avg_area_subset

upperlinfit = plot(t, ens_fit_sub[1, 1] .+ ens_fit_sub[2, 1] .* t, 
                    title = "Upper layer linear best fit of the growth\n of subset ensemble average area",
                    label = "Best fit of subset ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(t, ens_avg_area_subset[:, 1], 
    label = "Subset ensemble data")
#savefig(upperlinfit, "upperlinfitblob.png")

plot(t, ens_fit_sub[1, 2] .+ ens_fit_sub[2, 2] .* t, 
                    title = "Lower layer linear best fit of the growth\n of subset ensemble average area",
                    label = "Best fit of subset ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(t, ens_avg_area_subset[:, 2], 
    label = "Subset ensemble data")

K_sub_dim = @. ( ens_fit_sub[2, :] / (4 * π) ) * dims["Ld"] * 0.02

## Average area growth for members and ensemble average

member_subset_avg_area = first_moment(data, no_of_degrees * one_degree, 4)

first_mom_upper = plot(t, member_subset_avg_area[:, 1, 1], 
                        label = "Ensemble member", 
                        title = "(a) Upper layer average area growth for Gaussian blob initial condition,\ndata is subset every 60km meridionally and every "*string(no_of_degrees)*" degrees zonally",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
first_mom_lower = plot(t, member_subset_avg_area[:, 2, 1], 
                        label = "Ensemble member", 
                        title = "(b) Lower layer average area growth for Gaussian blob initial condition,\ndata is subset every 60km meridionally and every "*string(no_of_degrees)*" degrees zonally",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
for i ∈ 2:length(data)
    plot!(first_mom_upper, t, member_subset_avg_area[:, 1, i], label = false)
    plot!(first_mom_lower, t, member_subset_avg_area[:, 2, i], label = false)
end

plot!(first_mom_upper, t, ens_avg_area_subset[:, 1], label = "Ensemble average", line = (:dash, 2, :black))
plot!(first_mom_lower, t, ens_avg_area_subset[:, 2], label = "Ensemble average", line = (:dash, 2, :black))

fullplot = plot(first_mom_upper, first_mom_lower, layout = (2, 1), size= (800, 800))

## Diffusivity of the ensemble members in linear growth phase
j = 50
first_mom_upper = plot(t, member_subset_avg_area[:, 1, j], 
                        label = "Ensemble member", 
                        title = "Upper layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
first_mom_lower = plot(t, member_subset_avg_area[:, 2, j], 
                        label = "Ensemble member", 
                        title = "Lower layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)

member_subset_linfit = [[ones(length(t[round(Int64, 3*end / 4):end])) t[round(Int64, 3*end / 4):end]] \ member_subset_avg_area[round(Int64, 3*end / 4):end, :, k] for k ∈ 1:length(data)]

plot!(first_mom_upper, t[round(Int64, 3*end / 4):end], member_subset_linfit[j][1, 1] .+ member_subset_linfit[j][2, 1] .* t[round(Int64, 3*end / 4):end], label = "Linear fit")
plot!(first_mom_lower, t[round(Int64, 3*end / 4):end], member_subset_linfit[j][1, 2] .+ member_subset_linfit[j][2, 2] .* t[round(Int64, 3*end / 4):end], label = "Linear fit")
plot(first_mom_upper, first_mom_lower, layout = (2, 1), size = (800, 800))

## Extract the slope from the linear fits and find dimensional diffusivity
member_subset_linfit = [[member_subset_linfit[k][2, 1] for k ∈ 1:length(data)] [member_subset_linfit[k][2, 2] for k ∈ 1:length(data)]]'
K_subset_member_linfit = member_subset_linfit ./ (4 * π)
K_subset_member_linfit_dim = @. K_subset_member_linfit * dims["Ld"] * 0.02

μ_members = mean(K_subset_member_linfit_dim, dims = 2)
σ_members = std(K_subset_member_linfit_dim, dims = 2)

## Upper layer histogram of member diffusivities
upper_diff_hist_blob = fit(Histogram, K_subset_member_linfit_dim[1, :])
upper_diff_hist_blob_norm = normalize(upper_diff_hist_blob; mode = :probability)
upper_diff_hist_blob_plot = plot(upper_diff_hist_blob,
                                xlabel = "Diffusivity m²s⁻¹ ", 
                                ylabel = "Number of members",
                                title = "Histogram of ensemble members \n binned by diffusivity (upper layer)",
                                label = false, 
                                legend = :topright)
scatter!(upper_diff_hist_blob_plot, [K_sub_dim[1]], [0], label = "Ensemble average\ndiffusivity")
scatter!(upper_diff_hist_blob_plot, [findmin(K_subset_member_linfit_dim[1, :])[1]], [0], label = "Member with\nminimum diffisivity")
scatter!(upper_diff_hist_blob_plot, [findmax(K_subset_member_linfit_dim[1, :])[1]], [0], label = "Member with\nmaximum diffisivity")

## Lower layer histogram of member diffusivities
lower_diff_hist_blob = fit(Histogram, K_subset_member_linfit_dim[2, :], nbins = 12)
lower_diff_hist_blob_norm = normalize(lower_diff_hist_blob; mode = :probability)
lower_diff_hist_blob_plot = plot(lower_diff_hist_blob,
                            xlabel = "Diffusivity m²s⁻¹ ", 
                            ylabel = "Number of members",
                            title = "Histogram of ensemble members \n binned by diffusivity (lower layer)",
                            label = false, 
                            legend = :topright)
scatter!(lower_diff_hist_blob_plot, [K_sub_dim[2]], [0], label = "Ensemble average\ndiffusivity")
scatter!(lower_diff_hist_blob_plot, [findmin(K_subset_member_linfit_dim[2, :])[1]], [0], label = "Member with\nminimum diffisivity")
scatter!(lower_diff_hist_blob_plot, [findmax(K_subset_member_linfit_dim[2, :])[1]], [0], label = "Member with\nmaximum diffisivity")

## bootstrap ensemble average for varying subsets of data

N = 1000 
sample = 30
sample_vec = Array{Int64}(undef, sample)
sample_vals = 1:length(data)
ens_diff_subset = Array{Float64}(undef, N, 2)
#=
for n ∈ 1:N

    sample_vec = StatsBase.sample(sample_vals, sample, replace = false)
    data_sample = data[sample_vec]
    ens_conc_sample = ensemble_concentration(data_sample)
    ensemble_avg_sample = first_moment(ens_conc_sample, no_of_degrees * one_degree, 4)
    ens_fit_sample = [ones(length(t)) t] \ ensemble_avg_sample
    ens_diff_subset[n, :] =  ens_fit_sample[2, :]   

end

samples_diff_subset =  ens_diff_subset ./ (4 * π)   
dims = nondim2dim(data[1])
samples_diff_subset = @. samples_diff_subset * dims["Ld"] * 0.02

file = "bootstrap_blob_subset_4degreeszonal.jld2"
jldopen(file, "a+") do path
    path["bootstap"] = samples_diff_subset
end
=#

## Bootstrap for subset every 4 degrees zonally
bootstrap_4degrees = load("bootstrap_blob_subset_4degreeszonal.jld2")
samples_diff_subset = bootstrap_4degrees["bootstap"]

μ_samples = mean(samples_diff_subset, dims = 1)
σ_samples = std(samples_diff_subset, dims = 1)

## Upper layer
upper_bootstrap_hist = fit(Histogram, samples_diff_subset[:, 1])
upper_bootstrap_hist = normalize(upper_bootstrap_hist; mode = :probability)

#Upper layer
bootstrap_members_hist_upper = plot(upper_diff_hist_blob_norm, 
                                    xlabel = "Diffusivity m²s⁻¹", 
                                    ylabel = "Proportion of members",
                                    label = "Ensemble members",
                                    #title = "Diffusivity of each ensemble member and the bootstrapped \nsamples for the upper layer of the Gaussian blob",
                                    size = (800, 600),
                                    legend = :topleft)
#savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
plot!(bootstrap_members_hist_upper, upper_bootstrap_hist, label = "Bootstrapped samples")
#savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
scatter!(bootstrap_members_hist_upper, [μ_members[1]], [0], 
        markersize = 6,
        label = "Average diffusivity of\nensemble members")
#savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
scatter!(bootstrap_members_hist_upper, [μ_members[1] - σ_members[1], μ_members[1] + σ_members[1]], [0, 0], 
        markersize = 6,
        label = "± one standard deviation\nof ensemble members")
#savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper_mean_sd.png")
scatter!(bootstrap_members_hist_upper, [μ_samples[1]], [0],
        marker = :star,
        markersize = 6,
        label = "Average diffusivity of\nbootstrap samples")
#savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
scatter!(bootstrap_members_hist_upper, [μ_samples[1] - σ_samples[1], μ_samples[1] + σ_samples[1]], [0, 0], 
        marker = :star,
        markersize = 6,
        label = "± one standard deviation\nof bootstrap samples")
#savefig(bootstrap_members_hist_upper, "upper_blob_mem_boot.png")


## Lower layer
lower_bootstrap_hist = fit(Histogram, samples_diff_subset[:, 2])
lower_bootstrap_hist = normalize(lower_bootstrap_hist; mode = :probability)

bootstrap_members_hist_lower = plot(lower_diff_hist_blob_norm, 
                                    xlabel = "Diffusivity m²s⁻¹", 
                                    ylabel = "Proportion of members",
                                    label = "Ensemble members",
                                    #title = "Diffusivity of each ensemble member and the bootstrapped \nsamples for the lower layer of the Gaussian blob",
                                    size = (800, 600))
plot!(bootstrap_members_hist_lower, lower_bootstrap_hist, label = "Bootstrapped samples")
scatter!(bootstrap_members_hist_lower, [μ_members[2]], [0], 
        markersize = 6,
        label = "Average diffusivity of\nensemble members")
scatter!(bootstrap_members_hist_lower, [μ_members[2] - σ_members[2], μ_members[2] + σ_members[2]], [0, 0], 
        markersize = 6,
        label = "± one standard deviation\nof ensemble members")
scatter!(bootstrap_members_hist_lower, [μ_samples[2]], [0],
        marker = :star,
        markersize = 6,
        label = "Average diffusivity of\nbootstrap samples")
scatter!(bootstrap_members_hist_lower, [μ_samples[2] - σ_samples[2], μ_samples[2] + σ_samples[2]], [0, 0], 
        marker = :star,
        markersize = 6,
        label = "± one standard deviation\nof bootstrap samples")

##################################################################################################################
## Alternate (but maybe not as good method with the subset of data) of finding diffusivity using rise over run
Δt_mem = t[end] - t[round(Int64, 3*end / 4)]
ΔA_mem = member_subset_avg_area[end, :, :] .- member_subset_avg_area[round(Int64, 3*end / 4), :, :]

K_member_subset_ens = ΔA_mem ./ (4 * π * Δt_mem)

K_member_subset_dim = @. K_member_subset_ens * dims["Ld"] * 0.02
K_member_subset_dim[:, j]
##