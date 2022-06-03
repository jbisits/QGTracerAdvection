# These are plots for the tracer mixing paper. 

using CairoMakie, JLD2, Statistics, GLM

cd(joinpath(pwd(), "Paper plots"))
SimPath = joinpath("..", "Honours thesis/Experiment")

# If needed
module_path = "/Users/Joey/Documents/GitHub/QGTracerAdvection/Modules"
include(joinpath(module_path, "MeasureMixing.jl"))
using .MeasureMixing
################################################################################################
# Diffusion experiments
################################################################################################

diff_expt_path = joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 16, nx = 128, save_freq = 50, IC = GaussianBlob, Ensemble = false/SimulationData.jld2")
diff_expt_data = load(diff_expt_path)

diff_expt_plot = Figure(resolution = (1200, 1200), fontsize = 17)

titles = ["(a) Initial time" "(b) Initial time"; "(c) Final time" "(d) Final time"]
xlabs = ["Accumulated area", L"\hat{x}"]
xscales = [log10, identity]
ylabs = [L"Concentration ($\hat{C})", L"\hat{y}"]
ax = [Axis(diff_expt_plot[i, j],
            title = titles[i, j],
            xlabel = xlabs[j],
            xscale = xscales[j],
            ylabel = ylabs[j]) for i ∈ 1:2, j ∈ 1:2]

x, y = diff_expt_data["grid/x"], diff_expt_data["grid/y"]
for i ∈ 0:1

    hm = CairoMakie.heatmap!(ax[3 + i], x, y, diff_expt_data["snapshots/Concentration/"*string(i * 7000)],
                colormap = :deep)
    Colorbar(diff_expt_plot[1 + i, 3], hm, label = L"Concentration $(\hat{C})$") # color bar hidden
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

## Linear fit to diffusion data
t_diff = time_vec(diff_expt_data)
t_mat = [ones(length(t_diff)) t_diff]
first_moms_diff_expt = reshape(first_moment(diff_expt_data), :)
diff_model = lm(t_mat, first_moms_diff_expt)
int, slope = coef(diff_model)
r2(diff_model)

slope / (4*π)

# Plotting to see if fit is linear.
lines!(ax, t, int .+ t .* slope)
av_area_de
save("av_area_de.png", av_area_de)
################################################################################################
# Tracer concentration throughout experiments
################################################################################################
conc_data = load("SimulationData.jld2")
dims = nondim2dim(conc_data)
dims["Lx"] / 1000 # 7645km is zonal and meridional length
(dims["Δt"] * 18000) / 86400 # 1555 days, 
((dims["Δt"] * 18000) / 86400) / 365 # 4.2 years

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
tracer_plots = Figure(resolution = (1200, 1400), fontsize = 16)
plot_steps = 0:3000:15000
plot_steps_mat = reshape(plot_steps, (2, 3))
plot_times = round.(Int, [conc_data["snapshots/t/"*string(i)] for i ∈ plot_steps])
plot_times = reshape(plot_times, (3, 2))'
x̂ = conc_data["grid/x"]
ŷ = conc_data["grid/y"]
conc_plot_data = [abs.(conc_data["snapshots/Concentration/"*string(plot_steps_mat[j, i])][:, :, 1]) for j ∈ 1:2, i ∈ 1:3]
plot_letters = ["(a)" "(b)" "(c)"; "(d)" "(e)" "(f)"]

newline = "\n"
ax = [Axis(tracer_plots[i, j],
        xlabel = L"\hat{x} %$(newline)",
        ylabel = L"\hat{y}",
        title = L"%$(plot_letters[i, j]) \quad \hat{t} = %$(string(plot_times[i, j]))",
        aspect = 1
        ) for j ∈ 1:3, i ∈ 1:2]

for (i, axis) in enumerate(ax)

    plot_data = conc_plot_data[i]
    CairoMakie.heatmap!(axis, x̂, ŷ, plot_data, colormap = :deep)

end

for i ∈ 1:2, j ∈ 1:3

    plot_data = conc_plot_data[i, j]
    clims = (minimum(plot_data), maximum(plot_data))
    #cticks = round.(range(minimum(plot_data), maximum(plot_data), length = 3); sigdigits = 2)
    Colorbar(tracer_plots[i, j][2, 1], label = L"Concentration ($\hat{C}$)",
            limits = clims, vertical = false, 
            flipaxis = false, colormap = :deep, ticklabelrotation = 45.0)

end
tracer_plots
save("tracer_plots.png", tracer_plots)
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

area_inc = Figure(resolution = (400, 400))

titles = ["(a) Upper layer - square increase", "(b) Upper layer - meridional increase"]
ax = [Axis(area_inc[1, 1],
    title = "Upper layer",
    xlabel = L"\hat{t}",
    ylabel = L"%$(lang) \hat{A} \rangle",
)]

for files ∈ files_square

    data = load(files)
    Lx = round(Int, data["grid/Lx"])
    t = time_vec(data)
    first_mom_square = first_moment(data)
    lines!(ax[1], t, first_mom_square[:, 1],
        label = L"L\hat{x} = L\hat{y} = %$(string(Lx)) "#=*string(Lx)=#)

end
axislegend(ax[1], position = :lt)
area_inc
save("area_inc_sq.png", area_inc)

for files ∈ files_merid

    data = load(files)
    Ly = round(Int, data["grid/Ly"])
    t = time_vec(data)
    first_mom_merid = first_moment(data)
    lines!(ax[2], t, first_mom_merid[:, 1], 
    label = "Lŷ = "*string(Ly))

end

for i ∈ 1:2
    axislegend(ax[i], position = :lt)
end
axislegend(ax[1], position = :lt)
area_inc
save("area_inc_sq.png", area_inc)
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

μᵤ, μₗ = mean(member_diffs[:, 1]), mean(member_diffs[:, 2])
σᵤ_mem, σₗ_mem = std(member_diffs[:, 1]), std(member_diffs[:, 2])

# variability compared to ``true" diffusivity
100 * (μᵤ - σᵤ_mem) / ens_av_diffs[1]
100 * (μᵤ + σᵤ_mem) / ens_av_diffs[1]

100 * (μₗ - σₗ_mem) / ens_av_diffs[2]
100 * (μₗ + σₗ_mem) / ens_av_diffs[2]

# R² for linear fit to ensemble average
lm_data = [ones(length(t)) t]
linear_mod = lm(lm_data, ens_av_first_mom[:, 1])
r2(linear_mod)
lm_data = [ones(length(t)) t]
linear_mod = lm(lm_data, ens_av_first_mom[:, 2])
r2(linear_mod)

## First moment in time plots
first_moms_plot = Figure(resolution = (1000, 1000), fontsize = 17)

plot_time = 1:length(t)
short_plot_time = 1:findfirst(t .> 20)

#titles = ["(a) Upper layer" "(b) Lower layer"]
titles = [L"(a) Upper layer $\hat{t} = 0 - 90" L"(b) Upper layer $\hat{t} = 0 - 20";  L"(c) Lower layer $\hat{t} = 0 - 90" L"(d) Lower layer $\hat{t} = 0 - 20"]
ax = [Axis(first_moms_plot[i, j], 
            xlabel = L"\hat{t}",
            ylabel = L"\langle \hat{A} \rangle",
            title = titles[i, j], 
            aspect = 1) for j ∈ 1:2, i ∈ 1:2]

for i ∈ 1:length(member_first_moms[1, 1, :])
    if i == 1
        lines!(ax[1], t[plot_time], member_first_moms[plot_time, 1, i], color = :grey, label = "Ensemble member")
        lines!(ax[3], t[plot_time], member_first_moms[plot_time, 2, i], color = :grey, label = "Ensemble member")
        lines!(ax[2], t[short_plot_time], member_first_moms[short_plot_time, 1, i], color = :grey, label = "Ensemble member")
        lines!(ax[4], t[short_plot_time], member_first_moms[short_plot_time, 2, i], color = :grey, label = "Ensemble member")
    else
        lines!(ax[1], t[plot_time], member_first_moms[plot_time, 1, i], color = :grey)
        lines!(ax[3], t[plot_time], member_first_moms[plot_time, 2, i], color = :grey)
        lines!(ax[2], t[short_plot_time], member_first_moms[short_plot_time, 1, i], color = :grey)
        lines!(ax[4], t[short_plot_time], member_first_moms[short_plot_time, 2, i], color = :grey)
    end
end

lines!(ax[1], t[plot_time], ens_av_first_mom[plot_time, 1], label = "Ensemble mean")
lines!(ax[1], t[plot_time], ens_fit[1, 1] .+ t[plot_time] .* ens_fit[2, 1], 
        linestyle = :dash, 
        linewidth = 4,
        label = "Linear fit")
lines!(ax[2], t[short_plot_time], ens_av_first_mom[short_plot_time, 1], label = "Ensemble mean")
lines!(ax[2], t[short_plot_time], ens_fit[1, 1] .+ t[short_plot_time] .* ens_fit[2, 1], 
        linestyle = :dash, 
        linewidth = 4,
        label = "Linear fit")
lines!(ax[3], t[plot_time], ens_av_first_mom[plot_time, 2], label = "Ensemble mean")
lines!(ax[3], t[plot_time], ens_fit[1, 2] .+ t[plot_time] .* ens_fit[2, 2], 
        linestyle = :dash, 
        linewidth = 4,
        label = "Linear fit")
lines!(ax[4], t[short_plot_time], ens_av_first_mom[short_plot_time, 2], label = "Ensemble mean")
lines!(ax[4], t[short_plot_time], ens_fit[1, 2] .+ t[short_plot_time] .* ens_fit[2, 2], 
        linestyle = :dash, 
        linewidth = 4,
        label = "Linear fit")

axislegend(ax[1]; position = :lt)        
#Legend(first_moms_plot[3, :], ax[1])
#rowsize!(first_moms_plot.layout, 1, Relative(0.425))
#rowsize!(first_moms_plot.layout, 2, Relative(0.425))
first_moms_plot
save("first_moms.png", first_moms_plot)

################################################################################################
# Histograms of diffusivity
################################################################################################
## Member diffusivity
member_diffs = load("saved_data.jld2")["Diffusivity/member_diffs"]
bootstrap_samples = load("saved_data.jld2")["Bootstrap/diff_samples"]
bootstrap_samples_v2 = load("saved_data.jld2")["Bootstrap/diff_samples_v2"]
σᵤ, σₗ = std(bootstrap_samples_v2[:, 1]), std(bootstrap_samples_v2[:, 2])
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
    hist!(ax[i], bootstrap_samples_v2[:, i], normalization = :probability,
        label = "Bootstrapped diffusivity of ensemble average\nconcentration field")
    scatter!(ax[i], [mean(member_diffs[:, i])], [0], 
        label = "Mean diffusivity of ensemble members",
        color = :red)
    scatter!(ax[i], 
            [mean(member_diffs[:, i]) - std(member_diffs[:, i]), mean(member_diffs[:, i]) + std(member_diffs[:, i])], [0, 0], 
            label = "Mean diffusivity ± σ of ensemble members",
            color = :green)
    scatter!(ax[i], [ens_av_diffs[i]], [0],
            label = "Assumed \"true\" diffusvity",
            color = :magenta,
            marker = :diamond,
            markersize = 10)

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

# As percentage
upper_spatial_per = @. 100 * upper_spatial_rms_error / ens_av_diffs[1]
lower_spatial_per = @. 100 * lower_spatial_rms_error / ens_av_diffs[2]

upper_ts_per = @. 100 * upper_ts_rms_error / ens_av_diffs[1]
lower_ts_per = @. 100 * lower_ts_rms_error / ens_av_diffs[2]

# As decimal
upper_spatial_dec = @. upper_spatial_rms_error / ens_av_diffs[1]
lower_spatial_dec = @. lower_spatial_rms_error / ens_av_diffs[2]

upper_ts_dec = @. upper_ts_rms_error / ens_av_diffs[1]
lower_ts_dec = @. lower_ts_rms_error / ens_av_diffs[2]

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

100 * (upper_tempoal_rms_error[end] - upper_tempoal_rms_error[1]) / upper_tempoal_rms_error[1]
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

#####################################################################################
## Upper layer error plots 
#####################################################################################
# These are just upper layer error plots which will be used in the paper

## Version one, all spatial subsets no non linear colour scale
upper_err_plot = Figure(resolution = (600, 1200))

spat_RMS = upper_err_plot[1, 1]
temp_RMS = upper_err_plot[2, 1]
spat_temp_RMS = upper_err_plot[3, 1]

ax = [Axis(spat_RMS[1, 1], 
            xlabel = "Zonal distance between\ndata samples (km)",
            xticks = zonal_subset,
            xtickformat = xs -> [string(x .* 15) for x ∈ zonal_subset],
            xticklabelrotation = 45.0,
            ylabel = "Meridional distance between\ndata samples (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            xscale = log2, 
            yscale = log2,
            title = "(a) RMS error for spatial subsets\nof tracer concetration data", 
            aspect = 1)]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset, meridional_subset, log10.(upper_spatial_rms_error))
Colorbar(spat_RMS[1, 1][1, 2], upper_spatial,
        label = "log10(RMS) error of diffusivity from\nensemble members compared to 𝒦 (m²s⁻¹)")

ax = [Axis(temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "log10(RMS) error of diffusivity from\nensemble members compared to 𝒦 (m²s⁻¹)",
            xscale = log2,
            title = "(b) RMS error for temporal subsets\nof tracer concetration data",
            aspect = 1)]

upper_temp = lines!(ax[1], time_inc, log10.(upper_tempoal_rms_error))

ax = [Axis(spat_temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            xscale = log2,
            yscale = log2,
            title = "(c) RMS error for spatio-temporal\nsubsets of tracer concetration data",
            aspect = 1,
            alignmode = Inside())]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc, meridional_subset, log10.(upper_ts_rms_error))
Colorbar(spat_temp_RMS[1, 1][1, 2], upper_spatio_temp, 
        label = "log10(RMS) error of diffusivity from\nensemble members compared to 𝒦 (m²s⁻¹)")

upper_err_plot


## Version two, remove the largest spatial subset

upper_err_plot = Figure(resolution = (600, 1200))

spat_RMS = upper_err_plot[1, 1]
temp_RMS = upper_err_plot[2, 1]
spat_temp_RMS = upper_err_plot[3, 1]

ax = [Axis(spat_RMS[1, 1], 
            xlabel = "Zonal distance between\ndata samples (km)",
            xticks = zonal_subset,
            xtickformat = xs -> [string(x .* 15) for x ∈ zonal_subset[1:end-1]],
            xticklabelrotation = 45.0,
            ylabel = "Meridional distance between\ndata samples (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset[1:end-1]],
            xscale = log2, 
            yscale = log2,
            title = "(a) RMS error for spatial subsets\nof tracer concetration data", 
            aspect = 1)]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset[1:end-1], meridional_subset[1:end-1], upper_spatial_rms_error[1:end-1, 1:end-1])
Colorbar(spat_RMS[1, 1][1, 2], upper_spatial, 
        label = "RMS error of diffusivity from ensemble\nmembers compared to 𝒦 (m²s⁻¹)")

ax = [Axis(temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "RMS error of diffusivity from ensemble\nmembers compared to 𝒦 (m²s⁻¹)",
            xscale = log2,
            title = "(b) RMS error for temporal subsets\nof tracer concetration data",
            aspect = 1)]

upper_temp = lines!(ax[1], time_inc, upper_tempoal_rms_error)

ax = [Axis(spat_temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset[1:end-1]],
            xscale = log2,
            yscale = log2,
            title = "(c) RMS error for spatio-temporal subsets\nof tracer concetration data",
            aspect = 1,
            alignmode = Inside())]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc, meridional_subset[1:end-1], upper_ts_rms_error[:, 1:end-1])
Colorbar(spat_temp_RMS[1, 1][1, 2], upper_spatio_temp, 
        label = "RMS error of diffusivity from ensemble members compared to 𝒦 (m²s⁻¹)")

upper_err_plot


## Version three, remove the two largest spatial subsets

upper_err_plot = Figure(resolution = (600, 1200))

spat_RMS = upper_err_plot[1, 1]
temp_RMS = upper_err_plot[2, 1]
spat_temp_RMS = upper_err_plot[3, 1]

ax = [Axis(spat_RMS[1, 1], 
            xlabel = "Zonal distance between\ndata samples (km)",
            xticks = zonal_subset,
            xtickformat = xs -> [string(x .* 15) for x ∈ zonal_subset[1:end-2]],
            xticklabelrotation = 45.0,
            ylabel = "Meridional distance between\ndata samples (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset[1:end-2]],
            xscale = log2, 
            yscale = log2,
            title = "(a) RMS error for spatial subsets\nof tracer concetration data", 
            aspect = 1)]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset[1:end-2], meridional_subset[1:end-2], log10.(upper_spatial_rms_error[1:end-2, 1:end-2]))
Colorbar(spat_RMS[1, 1][1, 2], upper_spatial, #scale = log2,
        label = "log(RMS) error of diffusivity from ensemble\nmembers compared to 𝒦 (m²s⁻¹)")

ax = [Axis(temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "log(RMS) error of diffusivity from ensemble\nmembers compared to 𝒦 (m²s⁻¹)",
            xscale = log2,
            title = "(b) RMS error for temporal subsets\nof tracer concetration data",
            aspect = 1)]

upper_temp = lines!(ax[1], time_inc, log10.(upper_tempoal_rms_error))

ax = [Axis(spat_temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset[1:end-2]],
            xscale = log2,
            yscale = log2,
            title = "(c) RMS error for spatio-temporal subsets\nof tracer concetration data",
            aspect = 1,
            alignmode = Inside())]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc, meridional_subset[1:end-2], log10.(upper_ts_rms_error[:, 1:end-2]))
Colorbar(spat_temp_RMS[1, 1][1, 2], upper_spatio_temp, #scale = log10,
        label = "log(RMS) error of diffusivity from ensemble\nmembers compared to 𝒦 (m²s⁻¹)")

upper_err_plot
################################################################################################
# Using percentage of true diffusivity
################################################################################################
## Version one, all spatial subsets no non linear colour scale as percentage
upper_err_plot = Figure(resolution = (600, 1200))

spat_RMS = upper_err_plot[1, 1]
temp_RMS = upper_err_plot[2, 1]
spat_temp_RMS = upper_err_plot[3, 1]

ax = [Axis(spat_RMS[1, 1], 
            xlabel = "Zonal distance between\ndata samples (km)",
            xticks = zonal_subset,
            xtickformat = xs -> [string(x .* 15) for x ∈ zonal_subset],
            xticklabelrotation = 45.0,
            ylabel = "Meridional distance between\ndata samples (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            xscale = log2, 
            yscale = log2,
            title = "(a) RMS error for spatial subsets\nof tracer concetration data", 
            aspect = 1)]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset, meridional_subset, log10.(upper_spatial_per))
Colorbar(spat_RMS[1, 1][1, 2], upper_spatial,
        label = "log10(RMS) error of diffusivity\nas a percentage of 𝒦")

ax = [Axis(temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "log10(RMS) error of diffusivity\nas a percentage of 𝒦",
            xscale = log2,
            title = "(b) RMS error for temporal subsets\nof tracer concetration data",
            aspect = 1)]

upper_temp = lines!(ax[1], time_inc, log10.(100 .* upper_tempoal_rms_error ./ ens_av_diffs[1]))

ax = [Axis(spat_temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = time_inc,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = meridional_subset,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            xscale = log2,
            yscale = log2,
            title = "(c) RMS error for spatio-temporal\nsubsets of tracer concetration data",
            aspect = 1,
            alignmode = Inside())]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], time_inc, meridional_subset, log10.(upper_ts_per))
Colorbar(spat_temp_RMS[1, 1][1, 2], upper_spatio_temp,
        label = "log10(RMS) error of diffusivity\nas a percentage of 𝒦")

upper_err_plot
################################################################################################
# Using RMSe / 𝒦ᵤ as the scale
################################################################################################
## As decimal
upper_spatial_dec = @. upper_spatial_rms_error / ens_av_diffs[1]
upper_temporal_dec = @. upper_tempoal_rms_error ./ ens_av_diffs[1]
upper_ts_dec = @. upper_ts_rms_error / ens_av_diffs[1]

zonal_subset_nl = 0:8
meridional_subset_nl = 0:8
temporal_subset_nl = 0:6

upper_err_plot = Figure(resolution = (600, 1200))

spat_RMS = upper_err_plot[1, 1]
temp_RMS = upper_err_plot[2, 1]
spat_temp_RMS = upper_err_plot[3, 1]

ax = [Axis(spat_RMS[1, 1], 
            xlabel = "Zonal distance between\ndata samples (km)",
            xticks = zonal_subset_nl,
            xtickformat = xs -> [string(x .* 15) for x ∈ zonal_subset],
            xticklabelrotation = 45.0,
            ylabel = "Meridional distance between\ndata samples (km)",
            yticks = meridional_subset_nl,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            title = "(a) RMS error for spatial subsets\nof tracer concetration data", 
            aspect = 1)]

upper_spatial = CairoMakie.heatmap!(ax[1], zonal_subset_nl, meridional_subset_nl, log10.(upper_spatial_dec))

spatial_grid_vals = [Point(x, y) for x ∈ zonal_subset_nl for y ∈ meridional_subset_nl]
spatial_vals = string.(round.(reshape(100 * upper_spatial_dec', :); digits = 1))
text!(ax[1], spatial_vals, position = spatial_grid_vals, align = (:center, :center), color = :red, textsize = 10)

Colorbar(spat_RMS[1, 1][1, 2], upper_spatial,
        label = "log10(RMSe / 𝒦ᵤ)")

yidx = [1, 5, 6, 7]
ax = [Axis(temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = temporal_subset_nl,
            xtickformat = ts -> [string(t * 8) for t ∈ time_inc],
            yticks = round.(log10.(upper_temporal_dec)[yidx]; digits = 2),
            #ytickformat =  ts -> [string(round(100 *(10 ^ y); digits = 2)) for y ∈ temporal_ticks],
            ylabel = "log10(RMSe / 𝒦ᵤ)",
            title = "(b) RMS error for temporal subsets\nof tracer concetration data",
            aspect = 1)]

ax2 = [Axis(temp_RMS[1, 1], 
        #xlabel = "Time between data sampling (days)",
        #xticks = temporal_subset_nl,
        #xtickformat = ts -> [string(t * 8) for t ∈ time_inc],
        yticks = round.(log10.(upper_temporal_dec)[yidx]; digits = 2),
        ytickformat =  ts -> [string(round(100 * y; digits = 1)) for y ∈ upper_temporal_dec[yidx]],
        #ylabel = "(RMSe / 𝒦ᵤ) * 100",
        yticklabelcolor = :red,
        yaxisposition = :right,
        #title = "(b) RMS error for temporal subsets\nof tracer concetration data",
        aspect = 1)]        
   
hidespines!(ax2[1])
hidexdecorations!(ax2[1])
hideydecorations!(ax2[1], ticks = false)
upper_temp = lines!(ax[1], temporal_subset_nl, log10.(upper_temporal_dec))
lines!(ax2[1], temporal_subset_nl, log10.(upper_temporal_dec), overdraw = false)

ax = [Axis(spat_temp_RMS[1, 1], 
            xlabel = "Time between data sampling (days)",
            xticks = temporal_subset_nl,
            xtickformat = ts -> [string(t .* 8) for t ∈ time_inc],
            ylabel = "Zonal and meridional distance\nbetween data saplmes (km)",
            yticks = meridional_subset_nl,
            ytickformat = ys -> [string(y .* 15) for y ∈ meridional_subset],
            title = "(c) RMS error for spatio-temporal\nsubsets of tracer concetration data",
            aspect = 1,
            alignmode = Inside())]

upper_spatio_temp = CairoMakie.heatmap!(ax[1], temporal_subset_nl, meridional_subset_nl, log10.(upper_ts_dec))

st_grid_vals = [Point(x, y) for x ∈ temporal_subset_nl for y ∈ meridional_subset_nl]
st_vals = string.(round.(reshape(100 * upper_ts_dec', :); digits = 1))
text!(ax[1], st_vals, position = st_grid_vals, align = (:center, :center), color = :red, textsize = 10)

Colorbar(spat_temp_RMS[1, 1][1, 2], upper_spatio_temp,
        label = "log10(RMSe / 𝒦ᵤ)")

upper_err_plot

## Adding text values to the heatmap

spatial_grid_vals = [Point(x, y) for x ∈ zonal_subset for y ∈ meridional_subset]
spatial_vals = string.(round.(reshape(100 * upper_spatial_dec, :); digits = 2))
#scatter!(ax[1], spatial_grid_vals, marker = :circle)

text!(spat_RMS[1, 1], spatial_vals, position = spatial_grid_vals, justification = :center, color = (:black, 0.5))

for i ∈ 1:length(spatial_vals)

    t = text!(ax[1], spatial_vals[i],
              position = spatial_grid_vals[i],
              align = (:center, :center)
    )

end


###################################
## Figuring `text!`
###################################
f = Figure(resolution = (800, 800))

points = [Point(x, y) .* 200 for x in 1:3 for y in 1:3]
ax = Axis(f[1, 1])
scatter!(ax, points, marker = :circle, markersize = 10px)

symbols = (:left, :center, :right)

for ((justification, halign), point) in zip(Iterators.product(symbols, symbols), points)

    t = text!(ax, "a\nshort\nparagraph",
        color = (:black, 0.5),
        position = point,
        align = (halign, :center),
        justification = justification)

    #bb = boundingbox(t)
    #wireframe!(f, bb, color = (:red, 0.2))
end

for (p, al) in zip(points[3:3:end], (:left, :center, :right))
    text!(ax, "align :" * string(al), position = p .+ (0, 80),
        align = (:center, :baseline))
end

for (p, al) in zip(points[7:9], (:left, :center, :right))
    text!(ax, "justification\n:" * string(al), position = p .+ (80, 0),
        align = (:center, :top), rotation = pi/2)
end

f

## Problem occurs on the log2 scale.. Maybe I can re do the above plots not on log scale
f = Figure(resolution = (800, 800))

points = [Point(x, y) for x ∈ zonal_subset for y ∈ meridional_subset]
ax = Axis(f[1, 1], xscale = log2, yscale = log2)
scatter!(ax, points, marker = :circle, markersize = 10px)

text!(ax, "Works?", position = (log2(64), log2(8)), textsize = 5)

f