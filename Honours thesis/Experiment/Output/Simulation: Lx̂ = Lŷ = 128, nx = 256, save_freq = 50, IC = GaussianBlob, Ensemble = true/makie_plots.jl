# I needed log2 scale and it is not available in Plots. 
# These plots were definitely worth it. If I get on a roll may do the rest here.

using CairoMakie

################################################################################################
# Diffusion experiments
################################################################################################

diff_expt_path = joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 16, nx = 128, save_freq = 50, IC = GaussianBlob, Ensemble = false/SimulationData.jld2")
diff_expt_data = load(diff_expt_path)

diff_expt_plot = Figure(resolution = (1200, 1200))

titles = ["(a) Initial time" "(b) Initial time"; "(c) Final time" "(d) Final time"]
xlabs = ["x̂", "A"]
ylabs = ["ŷ", "Concentration"]
ylims_set = []
ax = [Axis(diff_expt_plot[i, j],
            title = titles[i, j],
            xlabel = xlabs[j],
            ylabel = ylabs[j],
            ylims = ylims_set,
            aspect = 1) for i ∈ 1:2, j ∈ 1:2]

x, y = diff_expt_data["grid/x"], diff_expt_data["grid/y"]
for i ∈ 0:1

    hm = CairoMakie.heatmap!(ax[1 + i], x, y, diff_expt_data["snapshots/Concentration/"*string(i * 7000)],
                colormap = :deep)
    Colorbar(diff_expt_plot[1 + i, 1][1, 2], hm) # color bar hidden
    colsize!(diff_expt_plot.layout, 1, Aspect(1, 1.0))
    lines!(ax[3 + i], 1:length(reshape(diff_expt_data["snapshots/Concentration/"*string(i * 7000)], :)), sort(reshape(diff_expt_data["snapshots/Concentration/"*string(i * 7000)], :), rev = true))
    CairoMakie.ylims!(ax[3 + i], high = maximum(diff_expt_data["snapshots/Concentration/"*string(0)]))

end

diff_expt_plot
################################################################################################
# Tracer experiment results and linear fits
################################################################################################
## Load in the data
first_mom_diff_data = load("av_area_diffs_linfits.jld2")
t = first_mom_diff_data["t"]
member_first_moms = first_mom_diff_data["member_first_moms"]
ens_av_first_mom = first_mom_diff_data["ens_av_first_mom"]
ens_fit = first_mom_diff_data["ens_fit"]
member_diffs = first_mom_diff_data["member_diffs"]

bootstrap_samples = load("bootstrap_blob.jld2")["bootstap"]

## First moment in time plots
first_moms_plot = Figure(resolution = (1000, 1000))

titles = ["(a) Upper Layer" "(b) Upper layer"; "(c) Lower Layer" "(d) Lower layer"]
ax = [Axis(first_moms_plot[i, j], 
            xlabel = "t (non-dimensional)",
            ylabel = "⟨A⟩ (non-dimensional)",
            title = titles[i, j], 
            aspect = 1) for i ∈ 1:2, j ∈ 1:2]

for i ∈ 1:length(member_first_moms[1, 1, :])
    if i == 1
        lines!(ax[1], t, member_first_moms[:, 1, i], color = :grey, label = "Ensemble member")
        lines!(ax[2], t, member_first_moms[:, 2, i], color = :grey, label = "Ensemble member")
    else
        lines!(ax[1], t, member_first_moms[:, 1, i], color = :grey)
        lines!(ax[2], t, member_first_moms[:, 2, i], color = :grey)
    end
end

lines!(ax[1], t, ens_av_first_mom[:, 1], linestyle = :dash, linewidth = 6, color = :black, label = "Ensemble average")
lines!(ax[2], t, ens_av_first_mom[:, 2], linestyle = :dash, linewidth = 6, color = :black, label = "Ensemble average")

lines!(ax[3], t, ens_av_first_mom[:, 1], label = "Ensemble average")
lines!(ax[3], t, ens_fit[1, 1] .+ t .* ens_fit[2, 1], 
        linestyle = :dash, 
        linewidth = 3,
        label = "Linear fit")
lines!(ax[4], t, ens_av_first_mom[:, 2], label = "Ensemble average")
lines!(ax[4], t, ens_fit[1, 2] .+ t .* ens_fit[2, 2], 
        linestyle = :dash, 
        linewidth = 3,
        label = "Linear fit")

for i ∈ 1:4
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
################################################################################################
# Subset data plots
################################################################################################
## Load in the data

subset_data = load("subset_data_for_plotting.jld2")

time_inc = subset_data["time_inc"]
zonal_subset = subset_data["zonal_subset"]
meridional_subset = subset_data["meridional_subset"]
spatial_subset = subset_data["spatial_subset"]
upper_spatial_rms_error = subset_data["upper_spatial_rms_error"]
lower_spatial_rms_error = subset_data["lower_spatial_rms_error"]
upper_tempoal_rms_error = subset_data["upper_tempoal_rms_error"]
lower_tempoal_rms_error = subset_data["lower_tempoal_rms_error"]
upper_ts_rms_error = subset_data["upper_spatiotemp_rms_error"]
lower_ts_rms_error = subset_data["lower_spatiotemp_rms_error"]

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

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset, meridional_subset, upper_spatial_rms_error)
lower_spatial = CairoMakie.heatmap!(ax[2], zonal_subset, meridional_subset, lower_spatial_rms_error)
Colorbar(spatial[1, 2], upper_spatial, label = "RMS error of ensemble\nmembers to 𝔎 (m²s⁻¹)")
Colorbar(spatial[2, 2], lower_spatial, label = "RMS error of ensemble\nmembers to 𝔎 (m²s⁻¹)")

colsize!(spatial.layout, 1, Aspect(1, 1.0))
spatial

save("spatial.png", spatial)

## Temporal subset of the data

temporal = Figure(resolution = (600, 800))

ax = [Axis(temporal[i, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 4) for t ∈ time_inc],
            ylabel = "RMS error of ensemble\nmembers to 𝔎 (m²s⁻¹)",
            xscale = log2,
            title = titles[i],
            aspect = 1) for i ∈ 1:2]

upper_temp = lines!(ax[1], time_inc, upper_tempoal_rms_error)
lower_temp = lines!(ax[2], time_inc, lower_tempoal_rms_error)

temporal

save("temporal.png", temporal)

## Spatio temporal subset of the data

spatio_temp = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(spatio_temp[i, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 4) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = spatial_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ spatial_subset],
            xscale = log2,
            yscale = log2,
            title = titles[i],
            aspect = 1) for i ∈ 1:2]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc, spatial_subset, upper_ts_rms_error)
lower_spatio_temp = CairoMakie.heatmap!(ax[2], time_inc, spatial_subset, lower_ts_rms_error)
Colorbar(spatio_temp[1, 2], upper_spatio_temp, label = "RMS error of ensemble\nmembers to 𝔎 (m²s⁻¹)")
Colorbar(spatio_temp[2, 2], lower_spatio_temp, label = "RMS error of ensemble\nmembers to 𝔎 (m²s⁻¹)")

colsize!(spatio_temp.layout, 1, Aspect(1, 1.0))
spatio_temp

save("spatio_temp.png", spatio_temp)