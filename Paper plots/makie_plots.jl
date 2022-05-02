# These are plots for the tracer mixing paper.
# 

using CairoMakie, JLD2, Statistics

cd(joinpath(pwd(), "Paper plots"))

################################################################################################
# Diffusion experiments
################################################################################################

diff_expt_path = joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 16, nx = 128, save_freq = 50, IC = GaussianBlob, Ensemble = false/SimulationData.jld2")
diff_expt_data = load(diff_expt_path)

diff_expt_plot = Figure(resolution = (1200, 1200))

titles = ["(a) Initial time" "(b) Initial time"; "(c) Final time" "(d) Final time"]
xlabs = ["𝒜", "x̂", ]
xscales = [log10, identity]
ylabs = ["Concentration", "ŷ"]
ax = [Axis(diff_expt_plot[i, j],
            title = titles[i, j],
            xlabel = xlabs[j],
            xscale = xscales[j],
            ylabel = ylabs[j],
            ylims = ylims_set) for i ∈ 1:2, j ∈ 1:2]

x, y = diff_expt_data["grid/x"], diff_expt_data["grid/y"]
for i ∈ 0:1

    hm = CairoMakie.heatmap!(ax[3 + i], x, y, diff_expt_data["snapshots/Concentration/"*string(i * 7000)],
                colormap = :deep)
    Colorbar(diff_expt_plot[1 + i, 3], hm, label = "Concentration") # color bar hidden
    lines!(ax[1 + i], 1:length(reshape(diff_expt_data["snapshots/Concentration/"*string(i * 7000)], :)), sort(reshape(diff_expt_data["snapshots/Concentration/"*string(i * 7000)], :), rev = true))
    CairoMakie.ylims!(ax[1 + i], high = maximum(diff_expt_data["snapshots/Concentration/"*string(0)]))

end
colsize!(diff_expt_plot.layout, 2, Aspect(1, 1))
diff_expt_plot
save("diff_expt_plot.png", diff_expt_plot)

## Average area during diffusion experiment
t = time_vec(diff_expt_data)
first_moms_diff_expt = first_moment(diff_expt_data)

av_area_de = Figure(resolution = (400, 400))
ax = Axis(av_area_de[1, 1],
        xlabel = "t̂",
        ylabel = "⟨Â⟩")
lines!(ax, t, first_moms_diff_expt[:, 1],
    label = "Growth of area\nof tracer patch")
axislegend(ax, position = :rb)
av_area_de
save("av_area_de.png", av_area_de)

################################################################################################
# Tracer concentration throughout experiments
################################################################################################
conc_data = load("SimulationData.jld2")

t = [conc_data["snapshots/t/"*string(i)] for i ∈ 0:100:18000]

x̂ = conc_data["grid/x"]
ŷ = conc_data["grid/y"]
conc = conc_data["snapshots/Concentration/0"][:, :, 1]

IC_conc = Figure(resolution = (500, 400))
ax = Axis(IC_conc[1, 1],
            xlabel = L"\hat{x}",
            ylabel = L"\hat{y}")

upper_IC = CairoMakie.heatmap!(ax, x̂, ŷ, conc, colormap = :deep)
Colorbar(IC_conc[1, 2], upper_IC, label = L"\hat{C}")
colsize!(IC_conc.layout, 1, Aspect(1, 1.0))

IC_conc
save("IC_conc.png", IC_conc)

# Evolution of tracer patch
tracer_plots = Figure(resolution = (1200, 1200))
plot_steps = 0:3000:15000
plot_times = round.(Int, [conc_data["snapshots/t/"*string(i)] for i ∈ plot_steps])
plot_times = reshape(plot_times, (3, 2))'
x̂ = conc_data["grid/x"]
ŷ = conc_data["grid/y"]
conc_plot_data = [abs.(conc_data["snapshots/Concentration/"*string(plot_steps_mat[i, j])][:, :, 1]) for i ∈ 1:2, j ∈ 1:3]

ax = [Axis(tracer_plots[i, j],
        xlabel = L"\hat{x}",
        ylabel = L"\hat{y}",
        title = L"\hat{t} = %$(string(plot_times[i, j]))",
        aspect = 1
        ) for j ∈ 1:3, i ∈ 1:2]
for (i, axis) in enumerate(ax)

    plot_data = conc_plot_data[i]
    clims = (minimum(plot_data), maximum(plot_data))
    heatmap!(axis, x̂, ŷ, plot_data, colormap = :deep)

end
tracer_plots
################################################################################################
# Growth of area
################################################################################################

## Load in the data
merid_area_inc = joinpath(SimPath, "Output/Meridional area increase blob")
files_merid = [joinpath(merid_area_inc, "SimulationData_64_64.jld2"), joinpath(merid_area_inc, "SimulationData_64_128.jld2"), 
        joinpath(merid_area_inc, "SimulationData_64_256.jld2")]

square_inc = joinpath(SimPath, "Output/Square area increase blob")
files_square = [joinpath(square_inc, "SimulationData_32.jld2"), joinpath(square_inc, "SimulationData_64.jld2"),
        joinpath(square_inc, "SimulationData_128.jld2"), joinpath(square_inc, "SimulationData_256.jld2")]

area_inc = Figure(resolution = (800, 800))

titles = ["(a) Upper layer" "(c) Upper layer"; "(b) Lower layer" "(d) Lower layer" ]
ax = [Axis(area_inc[i, j],
    title = titles[i, j],
    xlabel = "t̂",
    ylabel = "⟨Â⟩",
) for i ∈ 1:2, j ∈ 1:2]

for i ∈ 1:2

    for files ∈ files_square

        data = load(files)
        Lx = round(Int, data["grid/Lx"])
        t = time_vec(data)
        first_mom_square = first_moment(data)
        lines!(ax[i], t, first_mom_square[:, i],
            label = "Lx̂ = Lŷ = "*string(Lx))

    end

end

for i ∈ 3:4

    for files ∈ files_merid

        data = load(files)
        Ly = round(Int, data["grid/Ly"])
        t = time_vec(data)
        first_mom_square = first_moment(data)
        lines!(ax[i], t, first_mom_square[:, i - 2], 
        label = "Lŷ = "*string(Ly))

    end

end
for i ∈ 1:4
    axislegend(ax[i], position = :lt)
end
area_inc
save("area_inc.png", area_inc)
################################################################################################
# Tracer experiment results and linear fits
################################################################################################
## Load in the data
t = load("saved_data.jld2")["First_moms/t"]
member_first_moms = load("saved_data.jld2")["First_moms/Ensemble_first_moms"]
ens_av_first_mom = load("saved_data.jld2")["First_moms/Ensemble_avg_first_moms"]
ens_fit = load("saved_data.jld2")["First_moms/Ensemble_avg_lf"]
member_diffs = load("saved_data.jld2")["Diffusivity/member_diffs"]
ens_av_diffs = load("saved_data.jld2")["Diffusivity/ens_avg_diff" ]

bootstrap_samples = load("saved_data.jld2")["Bootstrap/diff_samples"]

## First moment in time plots
first_moms_plot = Figure(resolution = (1000, 1000))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(first_moms_plot[i, 1], 
            xlabel = "t̂",
            ylabel = "⟨Â⟩",
            title = titles[i], 
            aspect = 1) for i ∈ 1:2]

for i ∈ 1:length(member_first_moms[1, 1, :])
    if i == 1
        lines!(ax[1], t, member_first_moms[:, 1, i], color = :grey, label = "Ensemble member")
        lines!(ax[2], t, member_first_moms[:, 2, i], color = :grey, label = "Ensemble member")
    else
        lines!(ax[1], t, member_first_moms[:, 1, i], color = :grey)
        lines!(ax[2], t, member_first_moms[:, 2, i], color = :grey)
    end
end

lines!(ax[1], t, ens_av_first_mom[:, 1], label = "Ensemble average")
lines!(ax[1], t, ens_fit[1, 1] .+ t .* ens_fit[2, 1], 
        linestyle = :dash, 
        linewidth = 4,
        label = "Linear fit")
lines!(ax[2], t, ens_av_first_mom[:, 2], label = "Ensemble average")
lines!(ax[2], t, ens_fit[1, 2] .+ t .* ens_fit[2, 2], 
        linestyle = :dash, 
        linewidth = 4,
        label = "Linear fit")

for i ∈ 1:2
    axislegend(ax[i], position = :rb)
end

first_moms_plot
#save("first_moms.png", first_moms_plot)

################################################################################################
# Histograms of diffusivity
################################################################################################
## Member diffusivity
diffs_hist = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(diffs_hist[i, 1], 
            xlabel = "Diffusivity (m²s⁻¹)",
            ylabel = "Prortion of members",
            title = titles[i]) for i ∈ 1:2]

for i ∈ 1:2

    hist!(ax[i], member_diffs[:, i], normalization = :probability, bins = 10)
    CairoMakie.scatter!(ax[i], [mean(member_diffs[:, i])], [0], 
                    label = "Mean diffusivity of ensemble members",
                    color = :red)
    CairoMakie.scatter!(ax[i], 
                        [minimum(member_diffs[:, i])], [0], 
                        label = "Minimum diffusivity of ensemble members",
                        color = :orange)
    CairoMakie.scatter!(ax[i], 
                        [maximum(member_diffs[:, i])], [0], 
                        label = "Maximum diffusivity of ensemble members",
                        color = :green)

end
Legend(diffs_hist[3, 1], ax[1])
diffs_hist
save("diffs_hist.png", diffs_hist)
## Member diffusivity and bootstrap samples
bootstrap_hist = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(bootstrap_hist[i, 1], 
            xlabel = "Diffusivity (m²s⁻¹)",
            ylabel = "Prortion of members",
            title = titles[i]) for i ∈ 1:2]

for i ∈ 1:2

    hist!(ax[i], member_diffs[:, i], normalization = :probability, bins = 10,
        label = "Diffusivity of ensemble members")
    CairoMakie.scatter!(ax[i], [mean(member_diffs[:, i])], [0], 
                    label = "Mean diffusivity of ensemble members",
                    color = :red)
    CairoMakie.scatter!(ax[i], 
                        [mean(member_diffs[:, i]) - std(member_diffs[:, i]), mean(member_diffs[:, i]) + std(member_diffs[:, i])], [0, 0], 
                        label = "Mean diffusivity ± σ of ensemble members",
                        color = :green)
    hist!(ax[i], bootstrap_samples[:, i], normalization = :probability,
        label = "Bootstrapped diffusivity of ensemble average\nconcentration field")

end

Legend(bootstrap_hist[3, 1], ax[1])
bootstrap_hist
save("bootstrap_hist.png", bootstrap_hist)

## Percentage errors
((mean(member_diffs[:, 1]) - std(member_diffs[:, 1])), (mean(member_diffs[:, 1]) + std(member_diffs[:, 1])))
((mean(member_diffs[:, 2]) - std(member_diffs[:, 2])), (mean(member_diffs[:, 2]) + std(member_diffs[:, 2])))
upper_layer_per_err = 100 .* ((mean(member_diffs[:, 1]) - std(member_diffs[:, 1])) / ens_av_diffs[1], (mean(member_diffs[:, 1]) + std(member_diffs[:, 1])) / ens_av_diffs[1])
lower_layer_per_err = 100 .* ((mean(member_diffs[:, 2]) - std(member_diffs[:, 2])) / ens_av_diffs[2], (mean(member_diffs[:, 2]) + std(member_diffs[:, 2])) / ens_av_diffs[2])

################################################################################################
# Subset data plots
################################################################################################
## Load in the data

time_inc = 2 .^ (0:6)
zonal_subset = 2 .^ (0:8)
meridional_subset = 2 .^ (0:8)
spatial_subset = load("saved_data.jld2")["Spatial_subset/RMS_error"]
upper_spatial_rms_error = reshape(spatial_subset[:, 1, :], length(zonal_subset), length(meridional_subset))
lower_spatial_rms_error = reshape(spatial_subset[:, 2, :], length(zonal_subset), length(meridional_subset))
temporal_subset = load("saved_data.jld2")["Temporal_subset/RMS_error"]
upper_tempoal_rms_error = reshape(temporal_subset[:, 1, :], length(time_inc))
lower_tempoal_rms_error = reshape(temporal_subset[:, 2, :], length(time_inc))
spatio_temp_subset = load("saved_data.jld2")["Spatio_temp_subset/RMS_error"]
upper_ts_rms_error = reshape(spatio_temp_subset[:, 1, :, :], length(time_inc), length(meridional_subset))
lower_ts_rms_error = reshape(spatio_temp_subset[:, 2, :, :], length(time_inc), length(meridional_subset))

upper_spatial_per = @. 100 * upper_spatial_rms_error / ens_av_diffs[1]
lower_spatial_per = @. 100 * lower_spatial_rms_error / ens_av_diffs[2]

upper_ts_per = @. 100 * upper_ts_rms_error / ens_av_diffs[1]
lower_ts_per = @. 100 * lower_ts_rms_error / ens_av_diffs[2]
## Spatial subset of the data
spatial = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(spatial[i, 1], 
            xlabel = "Zonal distance between\ndata samples (km)",
            xticks = zonal_subset,
            xtickformat = xs -> [string(x .* 15) for x ∈ zonal_subset],
            xticklabelrotation = 45.0,
            ylabel = "Meridional distance between\ndata samples (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            xscale = log2, 
            yscale = log2,
            title = titles[i], 
            aspect = 1) for i ∈ 1:2]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset, meridional_subset, upper_spatial_per)
lower_spatial = CairoMakie.heatmap!(ax[2], zonal_subset, meridional_subset, lower_spatial_per)
Colorbar(spatial[1, 2], upper_spatial, label = "RMS error as percentage of 𝒦")
Colorbar(spatial[2, 2], lower_spatial, label = "RMS error as percentage of 𝒦")

colsize!(spatial.layout, 1, Aspect(1, 1.0))
spatial

save("spatial.png", spatial)

## Temporal subset of the data

temporal = Figure(resolution = (600, 800))

ax = [Axis(temporal[i, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "RMS error of diffusivity from ensemble\nmembers compared to 𝒦 (m²s⁻¹)",
            xscale = log2,
            title = titles[i],
            aspect = 1) for i ∈ 1:2]

upper_temp = lines!(ax[1], time_inc, upper_tempoal_rms_error)
lower_temp = lines!(ax[2], time_inc, lower_tempoal_rms_error)

temporal

save("temporal.png", temporal)

# RMS percentage increase

(upper_tempoal_rms_error[end] - upper_tempoal_rms_error[1]) / upper_tempoal_rms_error[1]
(lower_tempoal_rms_error[end] - lower_tempoal_rms_error[1]) / lower_tempoal_rms_error[1]

## Spatio temporal subset of the data

spatio_temp = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(spatio_temp[i, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = spatio_temp,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            xscale = log2,
            yscale = log2,
            title = titles[i],
            aspect = 1) for i ∈ 1:2]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc, meridional_subset, upper_ts_per)
lower_spatio_temp = CairoMakie.heatmap!(ax[2], time_inc, meridional_subset, lower_ts_per)
Colorbar(spatio_temp[1, 2], upper_spatio_temp, label = "RMS error as percentage of 𝒦")
Colorbar(spatio_temp[2, 2], lower_spatio_temp, label =  "RMS error as percentage of 𝒦")

colsize!(spatio_temp.layout, 1, Aspect(1, 1.0))
spatio_temp

save("spatio_temp.png", spatio_temp)

# Percentage increases

@. round(Int, 100 * (upper_ts_rms_error[end, :] - upper_ts_rms_error[1, :]) / upper_ts_rms_error[1, :])
@. round(Int, 100 * (lower_ts_rms_error[end, :] - lower_ts_rms_error[1, :]) / lower_ts_rms_error[1, :])

## Check reshape

A = Array{Float64}(undef, 10, 2, 3)
for i ∈ 1:3
    A[:, :, i] = [1:10 1:10]
end
A
reshape(A[:, 1, :], 10, 3)

spatial_subset