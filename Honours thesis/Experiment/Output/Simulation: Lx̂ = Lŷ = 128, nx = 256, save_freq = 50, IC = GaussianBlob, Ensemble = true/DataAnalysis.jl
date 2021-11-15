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
first_mom_upper = plot(t, first_moms[:, 1, 1], 
                        label = "Ensemble member", 
                        title = "Upper layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
first_mom_lower = plot(t, first_moms[:, 2, 1], 
                        label = "Ensemble member", 
                        title = "Lower layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
for i ∈ 2:length(data)
    plot!(first_mom_upper, t, first_moms[:, 1, i], label = false)
    plot!(first_mom_lower, t, first_moms[:, 2, i], label = false)
end

ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

plot!(first_mom_upper, t, ensemble_avg[:, 1], label = "Ensemble average", line = (:dash, 2, :black))
plot!(first_mom_lower, t, ensemble_avg[:, 2], label = "Ensemble average", line = (:dash, 2, :black))

fullplot = plot(first_mom_upper, first_mom_lower, layout = (2, 1), size= (800, 800))
#savefig(fullplot, "Gaussianblob_128dom_td6000.png")

## Diffusivity
Δt = t[41] - t[1]
ΔA = ensemble_avg[41, :] .- ensemble_avg[1, :]
K = ΔA ./ (4 * π * Δt)

## Linear fit
ens_fit = [ones(length(t)) t] \ ensemble_avg

upperlinfit = plot(t, ens_fit[1, 1] .+ ens_fit[2, 1] .* t, 
                    title = "Upper layer linear best fit of the growth of \n the ensemble average area",
                    label = "Best fit of ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(upperlinfit, t, ensemble_avg[:, 1], 
    label = "Ensemble data")
#savefig(upperlinfit, "upperlinfitblob.png")

lowerlinfit = plot(t, ens_fit[1, 2] .+ ens_fit[2, 2] .* t, 
                    title = "Lower layer linear best fit of the growth of \n the ensemble average area",
                    label = "Best fit of ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(lowerlinfit, t, ensemble_avg[:, 2], 
    label = "Ensemble data")
#savefig(lowerlinfit, "lowerlinfitblob.png")

K_linfit = ens_fit[2, :] ./ (4 * π)

dims = nondim2dim(data[1])
K_linfit_dim = @. K_linfit * dims["Ld"] * 0.02

## Tracer plots
tracer_plots = tracer_plot(data[1]; plot_freq = 500)
upperlayerblob = plot(tracer_plots[:, 1]..., size = (1400, 1400))
#savefig(upperlayerblob, "upperlayertracerblob.png")
final_blob = plot(tracer_plots[end, 1], size = (1000, 600))
savefig(final_blob, "final_blob.png")
lowerlayerblob = plot(tracer_plots[:, 2]..., size = (1400, 1400))
#savefig(lowerlayerblob, "lowerlayertracerblob.png")
upperlayerblobIC = plot(tracer_plots[1, 1], size = (800, 400))
#savefig(upperlayerblobIC, "blobIC.png")

## Transformed single plots 
final_sample = sort(abs.(reshape(data[1]["snapshots/Concentration/4000"][:, :, 1], :)), rev = true)
ordered_conc_data_sample = plot(area, final_sample, 
                        xlabel = "A", 
                        ylabel = "Concentration",
                        title = "Concentration data ordered higest to lowest",
                        label = false)
savefig(ordered_conc_data_sample, "ordered_conc_data_sample.png")

## Ensemble plots
ens_plots = tracer_plot(ensemble_conc; plot_freq = 500)
upperblobens = plot(ens_plots[:, 1]..., size = (1400, 1400))
last_plot = plot(ens_plots[end, 1], size = (1200, 500))
savefig(last_plot, "last_plot.png")
#savefig(upperblobens, "upperblobens.png")
lowerblobens = plot(ens_plots[:, 2]..., size = (1400, 1400))
#savefig(lowerblobens, "lowerblobens.png")

## Transfromed ensemble plots
initial_conc = sort(reshape(ensemble_conc["snapshots/Concentration/0"][:, :, 1], :), rev = true)
final_conc = sort(reshape(ensemble_conc["snapshots/Concentration/4000"][:, :, 1], :), rev = true)

area = 1:length(initial_conc)

plot(area, initial_conc)
ordered_conc_data = plot(area, final_conc, 
                        xlabel = "A", 
                        ylabel = "Concentration",
                        title = "Concentration data ordered higest to lowest",
                        label = false)


## Single and ensmeble plot for upper and lower layer
upper_ens_single = plot(tracer_plots[end, 1], ens_plots[end, 1], size = (1200, 600))
#savefig(upper_ens_single, "upper_ens_single.png")

## Diffusivity of each ensemble member 
j = 32
first_mom_upper = plot(t, first_moms[:, 1, j], 
                        label = "Ensemble member", 
                        title = "Upper layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
first_mom_lower = plot(t, first_moms[:, 2, j], 
                        label = "Ensemble member", 
                        title = "Lower layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)

Δt_mem = t[end] - t[round(Int64, 3*end / 4)]
ΔA_mem = first_moms[end, :, :] .- first_moms[round(Int64, 3*end / 4), :, :]

K_ens = ΔA_mem ./ (4 * π * Δt_mem)

dims = nondim2dim(data[1])
K_ens_dim = @. K_ens * dims["Ld"] * 0.02

#Upper layer
upper_diff_hist_blob = fit(Histogram, K_ens_dim[1, :], nbins = 12)
upper_diff_hist_blob_norm = normalize(upper_diff_hist_blob; mode = :probability)
upper_diff_hist_blob_plot = plot(upper_diff_hist_blob,
                                xlabel = "Diffusivity m²s⁻¹ ", 
                                ylabel = "Number of members",
                                title = "Histogram of ensemble members \n binned by diffusivity (upper layer)",
                                label = false, 
                                legend = :topright)
scatter!(upper_diff_hist_blob_plot, [K_linfit_dim[1]], [0], label = "Ensemble average\ndiffusivity")
scatter!(upper_diff_hist_blob_plot, [findmin(K_ens_dim[1, :])[1]], [0], label = "Member with\nminimum diffisivity")
scatter!(upper_diff_hist_blob_plot, [findmax(K_ens_dim[1, :])[1]], [0], label = "Member with\nmaximum diffisivity")
savefig(upper_diff_hist_blob_plot, "upper_diff_hist_blob.png")

#Lower layer
lower_diff_hist_blob = fit(Histogram, K_ens_dim[2, :], nbins = 8)
lower_diff_hist_blob_norm = normalize(lower_diff_hist_blob; mode = :probability)
lower_diff_hist_blob_plot = plot(lower_diff_hist_blob,
                            xlabel = "Diffusivity m²s⁻¹ ", 
                            ylabel = "Number of members",
                            title = "Histogram of ensemble members \n binned by diffusivity (upper layer)",
                            label = false, 
                            legend = :topright)
scatter!(lower_diff_hist_blob_plot, [K_linfit_dim[2]], [0], label = "Ensemble average\ndiffusivity")
scatter!(lower_diff_hist_blob_plot, [findmin(K_ens_dim[2, :])[1]], [0], label = "Member with\nminimum diffisivity")
scatter!(lower_diff_hist_blob_plot, [findmax(K_ens_dim[2, :])[1]], [0], label = "Member with\nmaximum diffisivity")
savefig(lower_diff_hist_blob_plot, "lower_diff_hist_blob.png")

#Summary stats 
μ_members = mean(K_ens_dim, dims = 2)
σ_members = std(K_ens_dim, dims = 2)

## Bootstrap the ensemble members to form error for ensemble average

N = 1000 
sample = 30
sample_vec = Array{Int64}(undef, sample)
sample_vals = 1:length(data)
ens_diff = Array{Float64}(undef, N, 2)
#= uncomment to run
for n ∈ 1:N

    sample_vec = StatsBase.sample(sample_vals, sample, replace = false)
    data_sample = data[sample_vec]
    ens_conc_sample = ensemble_concentration(data_sample)
    ensemble_avg_sample = first_moment(ens_conc_sample)
    ens_fit_sample = [ones(length(t)) t] \ ensemble_avg_sample
    ens_diff[n, :] =  ens_fit_sample[2, :]   

end
=#
samples_diff =  ens_diff ./ (4*π)   
dims = nondim2dim(data[1])
samples_diff = @. samples_diff * dims["Ld"] * 0.02

histogram(samples_diff[:, 1])
std(samples_diff[:, 1])
histogram(samples_diff[:, 2])
std(samples_diff[:, 2])


#This is a big thing to run so have saved to .jld2

file = "bootstrap_blob.jld2"
jldopen(file, "a+") do path
    path["bootstap"] = samples_diff
end

## Bootstrap results.
bootsrap = load("bootstrap_blob.jld2")
samples_diff = bootsrap["bootstap"]

μ_samples = mean(samples_diff, dims = 1)
σ_samples = std(samples_diff, dims = 1)

upper_bootstrap_hist = fit(Histogram, samples_diff[:, 1])
upper_bootstrap_hist = normalize(upper_bootstrap_hist; mode = :probability)

lower_bootstrap_hist = fit(Histogram, samples_diff[:, 2])
lower_bootstrap_hist = normalize(lower_bootstrap_hist; mode = :probability)

#Upper layer
bootstrap_members_hist_upper = plot(upper_diff_hist_blob_norm, 
                                    xlabel = "Diffusivity m²s⁻¹", 
                                    ylabel = "Proportion of members",
                                    label = "Ensemble members",
                                    #title = "Diffusivity of each ensemble member and the bootstrapped \nsamples for the upper layer of the Gaussian blob",
                                    size = (800, 600),
                                    legend = :topleft)
savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
plot!(bootstrap_members_hist_upper, upper_bootstrap_hist, label = "Bootstrapped samples")
savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
scatter!(bootstrap_members_hist_upper, [μ_members[1]], [0], 
        markersize = 6,
        label = "Average diffusivity of\nensemble members")
savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
scatter!(bootstrap_members_hist_upper, [μ_members[1] - σ_members[1], μ_members[1] + σ_members[1]], [0, 0], 
        markersize = 6,
        label = "± one standard deviation\nof ensemble members")
savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper_mean_sd.png")
scatter!(bootstrap_members_hist_upper, [μ_samples[1]], [0],
        marker = :star,
        markersize = 6,
        label = "Average diffusivity of\nbootstrap samples")
savefig(bootstrap_members_hist_upper, "bootstrap_members_hist_upper.png")
scatter!(bootstrap_members_hist_upper, [μ_samples[1] - σ_samples[1], μ_samples[1] + σ_samples[1]], [0, 0], 
        marker = :star,
        markersize = 6,
        label = "± one standard deviation\nof bootstrap samples")
savefig(bootstrap_members_hist_upper, "upper_blob_mem_boot.png")

n = 1:-0.1:0.1
upper_samples_inrange = Array{Int64}(undef, length(n))

for i ∈ 1:length(n)
    upper_samples_inrange[i] = length(findall(μ_members[1] -  n[i] * σ_members[1] .<= samples_diff[:, 1] .<= μ_members[1] + n[i] * σ_members[1]))
end
[n upper_samples_inrange]

# Percentage increase/decrease of ensemble average that captures all of the ensemble members.
# The upper layer captures all the samples from bootstrap in ± 0.7σ so we use this.

std_mul = 1
upper_diff = K_linfit_dim[1]
σ_upper = std_mul  * σ_members[1]
μ_upper = μ_members[1]
lower_lim, upper_lim = μ_upper - σ_upper, μ_upper + σ_upper
lower_per, upper_per = 100 * lower_lim / upper_diff, 100 *upper_lim / upper_diff

#Lower layer
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
savefig(bootstrap_members_hist_lower, "lower_blob_mem_boot.png")

n = 1:-0.1:0.1
lower_samples_inrange = Array{Int64}(undef, length(n))

for i ∈ 1:length(n)
    lower_samples_inrange[i] = length(findall(μ_members[2] -  n[i] * σ_members[2] .<= samples_diff[:, 2] .<= μ_members[2] + n[i] * σ_members[2]))
end
[n lower_samples_inrange]

# Percentage increase/decrease of ensemble average that captures all of the ensemble members.
# The upper layer captures all the samples from bootstrap in ± 0.7σ so we use this.

std_mul = 1
lower_diff = K_linfit_dim[2]
σ_lower = std_mul * σ_members[2]
μ_lower = μ_members[2]
lower_lim, upper_lim = μ_lower - σ_lower, μ_lower + σ_lower
lower_per, upper_per = 100 * lower_lim / lower_diff, 100 *upper_lim / lower_diff

##############################################################################################################
#Old/unused

## Or could fit Gaussian's and plot where the ensembel average diffusivity is
#Upper layer
upperlayer_normfit = fit(Normal, K_ens_dim[1, :])
upper_vals = minimum(K_ens_dim[1, :]):maximum(K_ens_dim[1, :])
upperlayer_normpdf = [pdf(upperlayer_normfit, x) for x ∈ upper_vals]

plot(upper_vals, upperlayer_normpdf, label = "Fitted normal to \nupper layer diffusivity\nestiamtes", legend = :topleft)
scatter!([K_linfit_dim[1]], [0], label = "Ensemble average\ndiffusivity")
scatter!([upperlayer_normfit.μ], [0], label = "Mean of diffusivity\nestimates")

#Lower layer
lowerlayer_normfit = fit(Normal, K_ens_dim[2, :])
lower_vals = minimum(K_ens_dim[2, :]):maximum(K_ens_dim[2, :])
lowerlayer_normpdf = [pdf(lowerlayer_normfit, x) for x ∈ lower_vals]

plot(lower_vals, lowerlayer_normpdf, label = "Fitted normal to \nlower layer diffusivity\nestiamtes", legend = :topright)
scatter!([K_linfit_dim[2]], [0], label = "Ensemble average\ndiffusivity")
scatter!([lowerlayer_normfit.μ], [0], label = "Mean of diffusivity\nestimates")


err_val = [100, 200, 300, 400, 500]
no_of_mems = Array{Float64}(undef, length(err_val))
for i in 1:length(err_val)
    lower_bound, upper_bound = K_linfit_dim[1] - err_val[i], K_linfit_dim[1] + err_val[i]
    temp = findall(K_ens_dim[1, :] .>= lower_bound)
    no_of_mems[i] = length(findall(K_ens_dim[1, temp] .<= upper_bound))
end

lower_bound, upper_bound = K_linfit_dim[1] - std(samples_diff[:, 1]),  K_linfit_dim[1] + std(samples_diff[:, 1])
temp = findall(K_ens_dim[1, :] .>= lower_bound)
findall(K_ens_dim[1, temp] .<= upper_bound)
findall(lower_bound .<= K_ens_dim[1, temp] .<= upper_bound)

err_val = σ_samples .* [1:10 1:10]
no_of_mems = similar(err_val)

for i in 1:length(err_val[:, 1])

    for j in 1:2
        lower_bound, upper_bound = K_linfit_dim[j] - err_val[i, j], K_linfit_dim[j] + err_val[i, j]
        no_of_mems[i, j] = length(findall(lower_bound .<=  K_ens_dim[j, :] .<= upper_bound))
    end
end

no_of_mems