# I needed log2 scale and it is not available in Plots. 
# These plots were definitely worth it. If I get on a roll may do the rest here.

# Note: No data is saved here. Need to compute it in the relevant scripts.
using CairoMakie

## Spatial subset of the data
spatial = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(spatial[i, 1], 
            xlabel = "Zonal subset",
            ylabel = "Meridional subset",
            xscale = log2, 
            yscale = log2,
            title = titles[i], 
            aspect = 1) for i ∈ 1:2]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset, meridional_subset, upper_av_err)
lower_spatial = CairoMakie.heatmap!(ax[2], zonal_subset, meridional_subset, lower_av_err)
Colorbar(spatial[1, 2], upper_spatial, label = "RMS error of ensemble\nmembers to κ")
Colorbar(spatial[2, 2], lower_spatial, label = "RMS error of ensemble\nmembers to κ")

colsize!(spatial.layout, 1, Aspect(1, 1.0))
spatial

## Spatio temporal subset of the data

spatio_temp = Figure(resolution = (600, 800))

titles = ["(a) Upper Layer", "(b) Lower layer"]
ax = [Axis(spatio_temp[i, 1], 
            xlabel = "Time between sample (days)",
            ylabel = "Distance between data samples",
            xscale = log2,
            yscale = log2,
            title = titles[i],
            aspect = 1) for i ∈ 1:2]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc .* 4, spatial_subset, upper_ts_rms_err)
lower_spatio_temp = CairoMakie.heatmap!(ax[2], time_inc .* 4, spatial_subset, lower_ts_rms_err)
Colorbar(spatio_temp[1, 2], upper_spatio_temp, label = "RMS error of ensemble\nmembers to κ")
Colorbar(spatio_temp[2, 2], lower_spatio_temp, label = "RMS error of ensemble\nmembers to κ")

colsize!(spatio_temp.layout, 1, Aspect(1, 1.0))
spatio_temp