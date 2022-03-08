# I needed log2 scale and it is not available in Plots. 
# These plots were definitely worth it. If I get on a roll may do the rest here.

using CairoMakie

## Load in the data

time_inc = load("subset_data_for_plotting.jld2")["time_inc"]
zonal_subset = load("subset_data_for_plotting.jld2")["zonal_subset"]
meridional_subset = load("subset_data_for_plotting.jld2")["meridional_subset"]
spatial_subset = load("subset_data_for_plotting.jld2")["spatial_subset"]
upper_spatial_rms_error = load("subset_data_for_plotting.jld2")["upper_spatial_rms_error"]
lower_spatial_rms_error = load("subset_data_for_plotting.jld2")["lower_spatial_rms_error"]
upper_tempoal_rms_error = load("subset_data_for_plotting.jld2")["upper_tempoal_rms_error"]
lower_tempoal_rms_error = load("subset_data_for_plotting.jld2")["lower_tempoal_rms_error"]
upper_ts_rms_error = load("subset_data_for_plotting.jld2")["upper_spatiotemp_rms_error"]
lower_ts_rms_error = load("subset_data_for_plotting.jld2")["lower_spatiotemp_rms_error"]

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