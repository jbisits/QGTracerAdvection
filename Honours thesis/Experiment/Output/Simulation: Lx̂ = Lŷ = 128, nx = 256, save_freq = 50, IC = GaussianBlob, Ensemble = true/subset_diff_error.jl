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
zonal_subset = 2 .^ (0:7)
meridional_subset = 2 .^ (0:7)
#=
# This is not that fast (though not that slow to run) but quicker to save and open. Of course can still change the way the data is subset.
member_diffs = Array{Float64}(undef, length(data), 2, length(zonal_subset) * length(meridional_subset))
k = 1
linear_time = 49 # This means there are 33 timesteps which are used for the time subsetting
linear_time_vec = t[linear_time:end]

for i ∈ zonal_subset, j ∈ meridional_subset
    
    first_moms = first_moment(data, i, j)
    member_subset_linfit = [[ones(length(linear_time_vec)) linear_time_vec] \ first_moms[linear_time:end, :, k] for k ∈ 1:length(data)]
    member_subset_linfit = [[member_subset_linfit[k][2, 1] for k ∈ 1:length(data)] [member_subset_linfit[k][2, 2] for k ∈ 1:length(data)]]
    K_subset_member_linfit = member_subset_linfit ./ (4π)
    member_diffs[:, :, k] = @. K_subset_member_linfit * dims["Ld"] * 0.02
    k += 1

end

file = "member_diffs_subset.jld2"
jldopen(file, "a+") do path
    path["member_diffs"] = member_diffs
end
=#
## Average absolute error heatmaps
member_diffs = load("member_diffs_subset.jld2")["member_diffs"]

## absolute error
member_diffs_err = @. abs(member_diffs - K_linfit_dim[1]) 
av_err = mean(member_diffs_err, dims = 1)

## RMS error
av_err = sqrt.( mean((member_diffs .- K_linfit_dim[1]).^2, dims = 1) )

upper_av_err = reshape(av_err[:, 1, :], (length(zonal_subset), length(meridional_subset)))
lower_av_err = reshape(av_err[:, 2, :], (length(zonal_subset), length(meridional_subset)))

zonal_points_per = 100 .* (zonal_subset ./ 256)
meridional_points_per = 100 .* (meridional_subset ./ 256)

upper_err_plot = heatmap(zonal_points_per, meridional_points_per, upper_av_err',
                    color = :viridis)

lower_err_plot = heatmap(zonal_points_per, meridional_points_per, lower_av_err',
                    color = :viridis)

err_plot = plot(upper_err_plot, lower_err_plot,
                    xlabel = "Percentage of zonal gridpoints",
                    xticks = [0, 5, 12, 25, 50],
                    ylabel = "Percentage of meridional gridpoints",
                    yticks = [0, 5, 12, 25, 50],
                    title = ["Upper layer absolute error for diffusivity\ncompared to ensemble average diffusivity" "Lower layer absolute error for diffusivity\ncompared to ensemble average diffusivity"],
                    colorbar_title = "RMS error of diffusivity (m²s⁻¹)",
                    color = :viridis,
                    layout = (2, 1), size = (1200, 1200),
                    #aspectratio = 1,
                    #framestyle = :box,
                    xlims = (0, zonal_points_per[end] + 12),
                    ylims = (0, meridional_points_per[end] + 12))

savefig(err_plot, "abs_error_heatmaps.png")

## Or can look at different subsets of the area on the heatmap
area_per = 100 .* ((zonal_subset .* meridional_subset) ./ 256^2)

upper_err_plot = heatmap(area_per, area_per, upper_av_err',
                    color = :viridis)

lower_err_plot = heatmap(area_per, area_per, lower_av_err',
                    color = :viridis)

err_plot = plot(upper_err_plot, lower_err_plot,
                    xlabel = "Percentage of zonal gridpoints",
                    xticks = ([0, 5, 12, 24], [0, 10, 25, 50]),
                    ylabel = "Percentage of meridional gridpoints",
                    yticks = ([0, 5, 12, 24], [0, 10, 25, 50]),
                    title = ["Upper layer absolute error for diffusivity\ncompared to ensemble average diffusivity" "Lower layer absolute error for diffusivity\ncompared to ensemble average diffusivity"],
                    colorbar_title = "RMS error of diffusivity (m²s⁻¹)",
                    color = :viridis,
                    layout = (2, 1), size = (1200, 1200),
                    aspectratio = 1,
                    framestyle = :box,
                    xlims = (0, area_per[end] + 7),
                    ylims = (0, area_per[end] + 7))

## Or can look at increase in square subsets only on a plot

upper_square_err = [upper_av_err[i, i] for i ∈ 1:length(upper_av_err[:, 1])]
lower_square_err = [lower_av_err[i, i] for i ∈ 1:length(lower_av_err[:, 1])]

area_per = 100 .* ((zonal_subset .* meridional_subset) ./ 256^2)

upper_err_square = plot(area_per, upper_square_err, label = false)
lower_err_square = plot(area_per, lower_square_err, label = false)
plot(upper_err_square, upper_err_square,
    xlabel = "Percentage of total area",
    ylabel = "RMS error from ensemble\nav. diff (m²s⁻¹)",
    title = ["Upper layer" "Lower layer"],
    layout = (2, 1),
    size = (1200, 1200))

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

k = 1
upper = plot(hist_diffs[k, 1], 
            xlabel = "Diffusivity m²s⁻¹",
            ylabel = "Proportion of members",
            title = "Upper Layer",
            label = "Diffusivity estimates")
scatter!(upper, [av_diff[k, 1]], [0], label = "Average of diffusivity estimates")
scatter!(upper, [av_diff[k, 1] + sd_diff[k, 1]], [0], label = "Av ± sd of diffusivity estimates", color = :green)
scatter!(upper, [av_diff[k, 1] - sd_diff[k, 1]], [0], label = false, color = :green)
plot!(upper, upper_bootstrap_hist, label = "Gold standard")

lower = plot(hist_diffs[k, 2], 
            xlabel = "Diffusivity m²s⁻¹",
            ylabel = "Proportion of members",
            title = "Upper Layer",
            label = "Diffusivity estimates")
scatter!(lower, [av_diff[k, 2]], [0], label = "Average of diffusivity estimates")
scatter!(lower, [av_diff[k, 2] + sd_diff[k, 2]], [0], label = "Av ± sd of diffusivity estimates", color = :green)
scatter!(lower, [av_diff[k, 2] - sd_diff[k, 2]], [0], label = false, color = :green)
plot!(lower, lower_bootstrap_hist, label = "Gold standard")

plot(upper, lower, size = (1000, 1000), layout = (2, 1))

## Now compute a diffusivity for all fifty ensemble members at each gridpoint in the subset of the domain.

square_subset = 2 .^ (0:1)
# This is not that fast (though not that slow to run) but quicker to save and open. Of course can still change the way the data is subset.
member_diffs = Array{Union{Missing, Float64}}(missing, length(data) * square_subset[end] * 2, 2, length(square_subset))
k = 1
linear_time = 49 # This means there are 33 timesteps which are used for the time subsetting
linear_time_vec = t[linear_time:end]

for i ∈ square_subset
    
    x_shift_vec = 0:(i - 1)
    y_shift_vec = 0:(i - 1)

    j = 0
    for m ∈ x_shift_vec, n ∈ y_shift_vec

        first_moms = first_moment(data, i, i, m, n)
        member_subset_linfit = [[ones(length(linear_time_vec)) linear_time_vec] \ first_moms[linear_time:end, :, k] for k ∈ 1:length(data)]
        member_subset_linfit = [[member_subset_linfit[k][2, 1] for k ∈ 1:length(data)] [member_subset_linfit[k][2, 2] for k ∈ 1:length(data)]]
        K_subset_member_linfit = member_subset_linfit ./ (4π)
        member_diffs[(1 + j * 50):(50 + j * 50), :, k] = @. K_subset_member_linfit * dims["Ld"] * 0.02
        j += 1
    
    end
    k += 1

end

member_diffs
member_diffs_rms_err = Array{Float64}(undef, 1, 2, length(member_diffs[1, 1, :]))

for i ∈ 1:length(member_diffs[1, 1, :])

    diffs = [collect(skipmissing(member_diffs[:, 1, i])) collect(skipmissing(member_diffs[:, 2, i]))]
    member_diffs_rms_err[:, :, i] = sqrt.( mean((diffs .- K_linfit_dim[1]).^ 2, dims = 1 ) )

end

member_diffs_rms_err

######### Varying temporal resolution of data
## First consider how it varies the full spatial data is used for differing teporal subsets
# Use 2ⁿ + 1 timesteps so always computing a beginning to end interval (if needed see my notes for explanation)
first_moms = first_moment(data)
linear_time = 49 # This means there are 33 timesteps which are used for the time subsetting

time_inc = @. 2^(0:5)
time_subset_linfits = Array{Float64}(undef, length(data), 2, length(time_inc))
k = 1
for i ∈ time_inc

    member_first_moms_linfit = [[ones(length(t[linear_time:i:end])) t[linear_time:i:end]] \ first_moms[linear_time:i:end, :, k] for k ∈ 1:length(data)]
    member_first_moms_linfit = [[member_first_moms_linfit[k][2, 1] for k ∈ 1:length(data)] [member_first_moms_linfit[k][2, 2] for k ∈ 1:length(data)]]
    time_subset_linfits[:, :, k] = member_first_moms_linfit ./ (4π)
    k += 1

end

dims = nondim2dim(data[1])
K_time_subset_linfits = @. time_subset_linfits * dims["Ld"] * 0.02

# Absolute error
#K_time_subset_linfits_rms_err = @. abs(K_time_subset_linfits - K_linfit_dim[1])
# RMS error
K_time_subset_linfits_rms_err = sqrt.( sum((K_time_subset_linfits .- K_linfit_dim[1]).^2, dims = 1) ./ length(K_time_subset_linfits[:, 1, 1]) )
#time_subset_av_err = mean(K_time_subset_linfits_rms_err, dims = 1)
upper_time_subset = reshape(K_time_subset_linfits_rms_err[:, 1, :], :)
lower_time_subset = reshape(K_time_subset_linfits_rms_err[:, 2, :], :)
plot(plot(time_inc .* 4, upper_time_subset, label = false), plot(time_inc .* 4, lower_time_subset, label = false), 
    xlabel = "Time increment (days)",
    ylabel = "RMS error of diffusivity (m²s⁻¹)",
    title = ["Upper layer" "Lower layer"],
    layout = (2, 1), 
    size = (800, 800))

## Histograms to see what is happening
j = 4
upper_hist = fit(Histogram, K_time_subset_linfits[:, 1, j], nbins = 12)
lower_hist = fit(Histogram, K_time_subset_linfits[:, 2, j], nbins = 12)
plot(plot(upper_hist), plot(lower_hist), 
    xlabel = "Diffusivity (m²s⁻¹)",
    ylabel = "Number of of members",
    title = ["Upper layer" "Lower layer"],
    layout = (2, 1), size = (800, 800))

######### Look at time varying but use more data 
## Use 2ⁿ steps to compute which allows more data but is not over the whole interval each time (again see notes for explanation)
# First just look at one observation from each time increment

first_moms = first_moment(data)
linear_time = 50 # This means there are 32 timesteps which are used for the time subsetting

time_inc = @. 2^(0:4)
time_subset_linfits = zeros(Float64, length(data), 2, length(time_inc))
k = 1

for i ∈ time_inc

    member_first_moms_linfit = [[ones(length(t[linear_time:i:end])) t[linear_time:i:end]] \ first_moms[linear_time:i:end, :, k] for k ∈ 1:length(data)]
    member_first_moms_linfit = [[member_first_moms_linfit[k][2, 1] for k ∈ 1:length(data)] [member_first_moms_linfit[k][2, 2] for k ∈ 1:length(data)]]
    time_subset_linfits[:, :, k] = member_first_moms_linfit ./ (4π)
    k += 1

end

dims = nondim2dim(data[1])
K_time_subset_linfits = @. time_subset_linfits * dims["Ld"] * 0.02

# RMS error
K_time_subset_linfits_rms_err = sqrt.( mean((K_time_subset_linfits .- K_linfit_dim[1]).^2, dims = 1) )
#time_subset_av_err = mean(K_time_subset_linfits_rms_err, dims = 1)
upper_time_subset = reshape(K_time_subset_linfits_rms_err[:, 1, :], :)
lower_time_subset = reshape(K_time_subset_linfits_rms_err[:, 2, :], :)
plot(plot(time_inc .* 4, upper_time_subset, label = false), plot(time_inc .* 4, lower_time_subset, label = false), 
    xlabel = "Time increment (days)",
    ylabel = "RMS error of diffusivity (m²s⁻¹)",
    title = ["Upper layer" "Lower layer"],
    layout = (2, 1), 
    size = (800, 800))

## Now taking diffusivity from all subsets within a time incremenet
# By including all different "starting points" of the vector can look at mulitple
# observations from each time subset.
first_moms = first_moment(data)
linear_time = 50 # This means there are 32 timesteps which are used for the time subsetting

time_inc = @. 2^(0:4)

time_subset_linfits = Array{Union{Missing, Float64}}(missing, time_inc[end] .* length(data), 2, length(time_inc))
k = 1
for i ∈ time_inc

    for j ∈ 0:(i - 1)

        member_first_moms_linfit = [[ones(length(t[(linear_time + j):i:end])) t[(linear_time + j):i:end]] \ first_moms[(linear_time + j):i:end, :, k] for k ∈ 1:length(data)]
        member_first_moms_linfit = [[member_first_moms_linfit[k][2, 1] for k ∈ 1:length(data)] [member_first_moms_linfit[k][2, 2] for k ∈ 1:length(data)]]
        time_subset_linfits[(1 + j * 50):(50 + j * 50), :, k] = member_first_moms_linfit ./ (4π)
    
    end

    k += 1

end

dims = nondim2dim(data[1])
K_time_subset_linfits = @. time_subset_linfits * dims["Ld"] * 0.02

# RMS error
K_time_subset_linfits_rms_err = Array{Float64}(undef, 1, 2, length(time_inc))

for i ∈ 1:length(time_inc)

    diffs = [collect(skipmissing(K_time_subset_linfits[:, 1, i])) collect(skipmissing(K_time_subset_linfits[:, 2, i]))]
    K_time_subset_linfits_rms_err[:, :, i] = sqrt.( mean((diffs .- K_linfit_dim[1]).^ 2, dims = 1 ) )

end

upper_time_subset = reshape(K_time_subset_linfits_rms_err[:, 1, :], :)
lower_time_subset = reshape(K_time_subset_linfits_rms_err[:, 2, :], :)
plot(plot(time_inc .* 4, upper_time_subset, label = false), plot(time_inc .* 4, lower_time_subset, label = false), 
    xlabel = "Time increment (days)",
    ylabel = "RMS error of diffusivity (m²s⁻¹)",
    title = ["Upper layer" "Lower layer"],
    layout = (2, 1), 
    size = (800, 800))

########### Time and spatial subsets
## Lastly look at subsets of spatial data that are then subset temporally    

linear_time = 50 # This means there are 32 timesteps which are used for the time subsetting
time_inc = @. 2^(0:4)
spatial_subset = @. 2^(0:7)

time_spatial_subset_linfits = Array{Union{Missing, Float64}}(missing, time_inc[end] .* length(data), 2, length(time_inc), length(spatial_subset))
n = 1
for m ∈ spatial_subset

    k = 1
    first_moms = first_moment(data, m, m)

    for i ∈ time_inc

        for j ∈ 0:(i - 1)

            member_first_moms_linfit = [[ones(length(t[(linear_time + j):i:end])) t[(linear_time + j):i:end]] \ first_moms[(linear_time + j):i:end, :, k] for k ∈ 1:length(data)]
            member_first_moms_linfit = [[member_first_moms_linfit[k][2, 1] for k ∈ 1:length(data)] [member_first_moms_linfit[k][2, 2] for k ∈ 1:length(data)]]
            time_spatial_subset_linfits[(1 + j * 50):(50 + j * 50), :, k, n] = member_first_moms_linfit ./ (4π)
        
        end

        k += 1

    end

    n += 1

end

time_spatial_subset_linfits_dims = @. time_spatial_subset_linfits * dims["Ld"] * 0.02
ts_rms_err = Array{Float64}(undef, 1, 2, length(time_inc), length(spatial_subset))

for j ∈ 1:length(spatial_subset)

    for i ∈ 1:length(time_inc)

        diffs = [collect(skipmissing(time_spatial_subset_linfits_dims[:, 1, i, j])) collect(skipmissing(time_spatial_subset_linfits_dims[:, 2, i, j]))]
        ts_rms_err[:, :, i, j] = sqrt.( mean((diffs .- K_linfit_dim[1]).^ 2, dims = 1 ) )

    end

end

ts_rms_err

upper_ts_rms_err = reshape(ts_rms_err[:, 1, :, :], (length(time_inc), length(spatial_subset)))
lower_ts_rms_err = reshape(ts_rms_err[:, 2, :, :], (length(time_inc), length(spatial_subset)))

area_per = @. 100 * (spatial_subset^2 / 256^2)

upper = heatmap(time_inc .* 4, area_per, upper_ts_rms_err', color = :viridis)
lower = heatmap(time_inc .* 4, area_per, lower_ts_rms_err', color = :viridis)

plot(upper, lower,
    xlabel = "Time subsample (days)",
    xticks = time_inc .* 4,
    ylabel = "Percentage of total area",
    yticks = [0, 1.5, 6.25, 25],
    colorbar_title = "RMS of concentration compared to ens. av (m²s⁻¹)",
    title = ["Upper layer" "Lower layer"],
    layout = (2, 1),
    size = (1200, 1200))