#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

## Load in the data for delay_time = Δt * 6000 with the new parameters
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
one_meridional_inc = 4 # 4 grid cells is one meridional increment again ≈ 60km

## Ensemble average diffusivity, our gold standard
ens_linfit = [ones(length(t)) t] \ ensemble_avg

K_linfit = ens_linfit[2, :] ./ (4 * π)

dims = nondim2dim(data[1])
K_linfit_dim = @. K_linfit * dims["Ld"] * 0.02


######### Varying spatial resolution of data

## Diffusivity estimates for ensemble members at different subsets of data
zonal_subset = 0:2:8
meridional_subset = 0:2:8
member_diffs = Array{Float64}(undef, length(data), 2, length(zonal_subset) * length(meridional_subset))
#= This is not that fast (though not that slow to run) but quicker to save and open. Of course can still change the way the data is subset.
k = 1
for i ∈ zonal_subset, j ∈ meridional_subset
    
    first_moms = first_moment(data, i * one_degree, j * one_meridional_inc)
    member_subset_linfit = [[ones(length(t[round(Int64, 3*end / 4):end])) t[round(Int64, 3*end / 4):end]] \ first_moms[round(Int64, 3*end / 4):end, :, k] for k ∈ 1:length(data)]
    member_subset_linfit = [[member_subset_linfit[k][2, 1] for k ∈ 1:length(data)] [member_subset_linfit[k][2, 2] for k ∈ 1:length(data)]]
    K_subset_member_linfit = member_subset_linfit ./ (4 * π)
    member_diffs[:, :, k] = @. K_subset_member_linfit * dims["Ld"] * 0.02
    k += 1

end

file = "member_diff_for_subsets.jld2"
jldopen(file, "a+") do path
    path["member_diffs"] = member_diffs
end
=#
## Average absolute error heatmaps

member_diffs = load("member_diff_for_subsets.jld2")["member_diffs"]
member_diffs_abs_err = @. abs(member_diffs - K_linfit_dim[1])

av_abs_err = mean(member_diffs_abs_err, dims = 1)

upper_av_err = reshape(av_abs_err[:, 1, :], (length(zonal_subset), length(meridional_subset)))
lower_av_err = reshape(av_abs_err[:, 2, :], (length(zonal_subset), length(meridional_subset)))

upper_err_plot = heatmap(zonal_subset .* 60, meridional_subset, upper_av_err', 
                    xlabel = "Meridional subset (km)",
                    ylabel = "Zonal subset (degrees)",
                    title = "Upper layer absolute error for diffusivity\ncompared to ensemble average diffusivity",
                    colorbar_title = "Average absolute error (m²s⁻¹)",
                    color = :viridis)

lower_err_plot = heatmap(zonal_subset .* 60, meridional_subset, lower_av_err', 
                    xlabel = "Meridional subset (km)",
                    ylabel = "Zonal subset (degrees)",
                    title = "Lower layer absolute error for diffusivity\ncompared to ensemble average diffusivity",
                    colorbar_title = "Average absolute error (m²s⁻¹)",
                    color = :viridis)

err_plot = plot(upper_err_plot, lower_err_plot, layout = (2, 1), size = (1200, 1200))

savefig(err_plot, "abs_error_heatmaps.png")

## Histograms of diffusivity estimates with coarser resolution and gold standard ensemble average diffusivity

bootsrap = load("bootstrap_blob.jld2")
samples_diff = bootsrap["bootstap"]

μ_samples = mean(samples_diff, dims = 1)
σ_samples = std(samples_diff, dims = 1)

upper_bootstrap_hist = fit(Histogram, samples_diff[:, 1])
upper_bootstrap_hist = normalize(upper_bootstrap_hist; mode = :probability)

lower_bootstrap_hist = fit(Histogram, samples_diff[:, 2])
lower_bootstrap_hist = normalize(lower_bootstrap_hist; mode = :probability)

# Member histograms for the different data subsets
# For now consider all the histograms and choose some (if any to show later)

hist_diffs = Array{Histogram}(undef, length(member_diffs[1, 1, :]), 2)
av_diff = Array{Float64}(undef, length(member_diffs[1, 1, :]), 2)
sd_diff = Array{Float64}(undef, length(member_diffs[1, 1, :]), 2)

for i ∈ 1:length(member_diffs[1, 1, :]), j ∈ 1:2
    
    hist = fit(Histogram, member_diffs[:, j, i])
    hist_diffs[i, j] = normalize(hist; mode = :probability)

    av_diff[i, j] = mean(member_diffs[:, j, i])
    sd_diff[i, j] = std(member_diffs[:, j, i])
end

k = 20
upper = plot(hist_diffs[k, 1], 
            xlabel = "Diffusivity m²s⁻¹",
            ylabel = "Proportion of members",
            title = "Upper Layer",
            label = "Diffusivity estimates")
scatter!(upper, [av_diff[k, 1]], [0], label = "Average of diffusivity estimates")
scatter!(upper, [av_diff[k, 1] + sd_diff[k, 1]], [0], label = "Av + sd of diffusivity estimates", color = :green)
scatter!(upper, [av_diff[k, 1] - sd_diff[k, 1]], [0], label = "Av - sd of diffusivity estimates", color = :green)
plot!(upper, upper_bootstrap_hist, label = "Gold standard")

lower = plot(hist_diffs[k, 2], 
            xlabel = "Diffusivity m²s⁻¹",
            ylabel = "Proportion of members",
            title = "Upper Layer",
            label = "Diffusivity estimates")
scatter!(lower, [av_diff[k, 2]], [0], label = "Average of diffusivity estimates")
scatter!(lower, [av_diff[k, 2] + sd_diff[k, 2]], [0], label = "Av + sd of diffusivity estimates", color = :green)
scatter!(lower, [av_diff[k, 2] - sd_diff[k, 2]], [0], label = "Av - sd of diffusivity estimates", color = :green)
plot!(lower, lower_bootstrap_hist, label = "Gold standard")

plot(upper, lower, size = (1000, 1000), layout = (2, 1))

######### Varying temporal resolution of data

first_moms = first_moment(data, 4 * one_degree, 4 * one_meridional_inc)
time_inc = 16
no_of_members = 1
plot(t[1:time_inc:end], [first_moms[1:time_inc:end, 1, i] for i ∈ 1:no_of_members], label = false)